# ===== Engineering Mathematics -- Question 21 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 2nd May, 2026

# Optimal Control: Disturbance Impact Comparison

using JuMP
using Ipopt
using Plots
using LinearAlgebra

g = 1.62

# ---- Optimization Function ---- #
function optimize_trajectory(
    N,
    state_zero,
    state_end,
    tspan,
    u_max,
    disturbance = false,
    method=:euler,
)   
    if (method != :euler)
        return 1
    end

    t_start, t_end = tspan
    h = (t_end - t_start) / N

    model = Model(Ipopt.Optimizer)
    set_silent(model)

    @variables(model, begin
        x[1:N+1]
        z[1:N+1]

        vx[1:N+1]
        vz[1:N+1]

        ux[1:N]  # Size N to match Forward Euler step intervals cleanly
        uz[1:N]
    end)

    # -- Solver Constraints -- #
    for k in 1:N
        if disturbance
            # Fixed profile sequence for this specific run instance
            α = rand() * (0.9 - 0.7) + 0.7
        else
            α = 1.0
        end

        @constraint(model, x[k+1] == x[k] + h * vx[k])
        @constraint(model, z[k+1] == z[k] + h * vz[k])

        @constraint(model, vx[k+1] == vx[k] + h * ux[k] * α)
        @constraint(model, vz[k+1] == vz[k] + h * (uz[k] * α - g))
    end

    # -- Define Objective -- #
    @objective(
        model,
        Min,
        sum(
            sqrt(ux[k]^2 + uz[k]^2 + 1e-6) * h for k in 1:N
        )
    )

    # -- Define Model Constraints -- #
    @constraint(model, x[1] == state_zero[1])
    @constraint(model, z[1] == state_zero[2])
    @constraint(model, vx[1] == state_zero[3])
    @constraint(model, vz[1] == state_zero[4])

    @constraint(model, x[N+1] == state_end[1])
    @constraint(model, z[N+1] == state_end[2])
    @constraint(model, vx[N+1] == state_end[3])
    @constraint(model, vz[N+1] == state_end[4])

    @constraint(model, z .>= 0)
    
    # Syntax cleaned: standard <= comparison operator and squared u_max bound
    @constraint(model, [k=1:N], ux[k]^2 + uz[k]^2 <= u_max^2)

    optimize!(model)

    time_array = range(t_start, t_end, length = N+1)
    return time_array, value.(x), value.(z), value.(vx), value.(vz), value.(ux), value.(uz)
end

# ---- Running Both Configurations ---- #
N = 100
state_zero = 5, 10, -1, -2
state_end = 0, 0, 0, 0
tspan = (0.0, 5.0)
u_max = 10.0

println("--- Running Nominal Case (No Disturbance) ---")
t, x_nom, z_nom, _, _, _, _ = optimize_trajectory(N, state_zero, state_end, tspan, u_max, false, :euler)

println("--- Running Perturbed Case (With Disturbance) ---")
_, x_dist, z_dist, _, _, _, _ = optimize_trajectory(N, state_zero, state_end, tspan, u_max, true, :euler)

# ---- Error Metrics and Computations ---- #

# Coordinate differences at each discrete time step
diff_x = abs.(x_dist .- x_nom)
diff_z = abs.(z_dist .- z_nom)

# Distance error metric at each step: Euclidean spatial separation 
step_distance_error = sqrt.(diff_x.^2 .+ diff_z.^2)

# Cumulative distance difference summed step-by-step
cumulative_difference = cumsum(step_distance_error)

# Final Position Assessment
final_err_x = abs(x_dist[end] - x_nom[end])
final_err_z = abs(z_dist[end] - z_nom[end])
final_euclidean_distance = sqrt(final_err_x^2 + final_err_z^2)

println("\n================ Boundary Results ================")
println("Nominal Final Position:   (", round(x_nom[end], digits=6), ", ", round(z_nom[end], digits=6), ")")
println("Perturbed Final Position: (", round(x_dist[end], digits=6), ", ", round(z_dist[end], digits=6), ")")
println("--------------------------------------------------")
println("Final X Coordinate Delta:  ", final_err_x, " m")
println("Final Z Coordinate Delta:  ", final_err_z, " m")
println("Total Final Position Gap:  ", final_euclidean_distance, " m (Solver Tolerance Convergence)")
println("==================================================")

# ---- Visualization Assembly ---- #

# Plot 1: Trajectory Profile Space Comparison (X vs Z)
p1 = plot(x_nom, z_nom, label="Nominal", lw=2.5, color=:blue, title="Spatial Trajectory Path (X vs Z)", xlabel="X [m]", ylabel="Z [m]")
plot!(p1, x_dist, z_dist, label="Perturbed", lw=2, linestyle=:dash, color=:red)

# Plot 2: Lateral Evolution Over Time (X vs t)
p2 = plot(t, x_nom, label="Nominal", lw=2, color=:blue, title="Lateral Position Profile", xlabel="t [s]", ylabel="X [m]")
plot!(p2, t, x_dist, label="Perturbed", lw=2, linestyle=:dash, color=:red)

# Plot 3: Altitude Evolution Over Time (Z vs t)
p3 = plot(t, z_nom, label="Nominal", lw=2, color=:blue, title="Altitude Profile", xlabel="t [s]", ylabel="Z [m]")
plot!(p3, t, z_dist, label="Perturbed", lw=2, linestyle=:dash, color=:red)

# Plot 4: Metric Tracking Graph (Cumulative Tracking Deviation vs Time)
p4 = plot(t, cumulative_difference, label="Σ ΔPosition", lw=2.5, color=:purple, legend=:topleft,
          title="Integrated Cumulative Path Deviation", xlabel="t [s]", ylabel="Accumulated Error [m]")

# Consolidate views into a clear 2x2 presentation grid
comparison_layout = plot(p1, p2, p3, p4, layout = (2, 2), size = (1100, 850))
display(comparison_layout)
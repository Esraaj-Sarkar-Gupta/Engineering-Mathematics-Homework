# ===== Engineering Mathematics -- Question 21 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 2nd May, 2026
# Optimal Control: Euler vs Trapezoid Comparison (Monte Carlo + Contour)

using JuMP
using Ipopt
using Plots
using Distributions
using LinearAlgebra

g = 1.62

println("---- Run ---- ")

# ---- Optimization Function ---- #
function optimize_trajectory(
    N,
    state_zero,
    state_end,
    tspan,
    u_max,
    method=:euler
)
    t_start, t_end = tspan
    h = (t_end - t_start) / N

    model = Model(Ipopt.Optimizer)
    set_silent(model)

    # -- Solver Constraints -- #
    if (method == :euler)
        @variables(model, begin
            Xx[1:N+1]; Xz[1:N+1]
            vx[1:N+1]; vz[1:N+1]
            ux[1:N]; uz[1:N]
        end)
        
        for k in 1:N
            @constraint(model, Xx[k+1] == Xx[k] + h * vx[k])
            @constraint(model, Xz[k+1] == Xz[k] + h * vz[k])
            @constraint(model, vx[k+1] == vx[k] + h * ux[k])
            @constraint(model, vz[k+1] == vz[k] + h * (uz[k] - g))
            # Energy/Control Constraint
            @constraint(model, ux[k]^2 + uz[k]^2 <= u_max^2)
        end

    elseif (method == :trapezoid)
        @variables(model, begin
            Xx[1:N+1]; Xz[1:N+1]
            vx[1:N+1]; vz[1:N+1]
            ux[1:N+1]; uz[1:N+1]
        end)

        for k in 1:N
            @constraint(model, Xx[k+1] == Xx[k] + h/2 * (vx[k] + vx[k+1]))
            @constraint(model, Xz[k+1] == Xz[k] + h/2 * (vz[k] + vz[k+1]))
            @constraint(model, vx[k+1] == vx[k] + h/2 * (ux[k] + ux[k+1]))
            @constraint(model, vz[k+1] == vz[k] + h/2 * (uz[k] + uz[k+1] - 2*g))
            # Energy/Control Constraint
            @constraint(model, ux[k]^2 + uz[k]^2 <= u_max^2)
        end
    else
        println("Optimization method was not recognized")
        return false, nothing
    end

    @objective(
        model,
        Min,
        sum(sqrt(ux[k]^2 + uz[k]^2 + 1e-6) * h for k in 1:N)
    )

    # Define Model Constraints
    @constraint(model, [Xx[1], Xz[1], vx[1], vz[1]] .== [state_zero...])
    @constraint(model, [Xx[N+1], Xz[N+1], vx[N+1], vz[N+1]] .== [state_end...])
    @constraint(model, Xz .>= 0)

    optimize!(model)
    status = termination_status(model)

    is_success = (status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL)
    return is_success
end


# ---- Monte Carlo Engine ---- #
function MC(num_samples, x_bounds, z_bounds)
    x_dist = Uniform(x_bounds[1], x_bounds[2])
    z_dist = Uniform(z_bounds[1], z_bounds[2])

    # Store raw MC data
    mc_x = Float64[]
    mc_z = Float64[]
    mc_status = Float64[] # 1.0 for feasible, 0.0 for infeasible

    println("Starting Monte Carlo Sampling ($num_samples points)...")

    for i in 1:num_samples
        sample_x = rand(x_dist)
        sample_z = rand(z_dist)
        
        state_zero = (sample_x, sample_z, -1.0, -2.0)
        state_end = (0.0, 0.0, 0.0, 0.0)

        # Using N=80 for slightly faster MC iterations
        is_feasible = optimize_trajectory(80, state_zero, state_end, (0.0, 5.0), 15.0, :trapezoid)

        push!(mc_x, sample_x)
        push!(mc_z, sample_z)
        push!(mc_status, is_feasible ? 1.0 : 0.0)
        
        if i % 50 == 0
            println("Completed $i / $num_samples samples...")
        end
    end

    println("Finished sampling! Interpolating data for contour plot...")
    return mc_x, mc_z, mc_status
end

# ---- Inverse Distance Weighting Interpolation ---- #
# Interpolates scattered MC data onto a regular grid for the contour plotter
function interpolate_to_grid(mc_x, mc_z, mc_status, bounds_x, bounds_z, resolution=50)
    grid_x = range(bounds_x[1], bounds_x[2], length=resolution)
    grid_z = range(bounds_z[1], bounds_z[2], length=resolution)
    grid_vals = zeros(resolution, resolution)

    for (i, gx) in enumerate(grid_x)
        for (j, gz) in enumerate(grid_z)
            weight_sum = 0.0
            val_sum = 0.0
            
            for k in 1:length(mc_x)
                dist = sqrt((gx - mc_x[k])^2 + (gz - mc_z[k])^2) + 1e-5 # <-- Avoid division by zero errors
                w = 1.0 / (dist^3) # Cubic weighting
                
                val_sum += w * mc_status[k]
                weight_sum += w
            end
            
            grid_vals[j, i] = val_sum / weight_sum
        end
    end
    
    return grid_x, grid_z, grid_vals
end

# ---- Execution ---- #
bounds_x = (0.0, 20.0)
bounds_z = (0.0, 20.0)
samples = 500

# Run Monte Carlo
mc_x, mc_z, mc_status = MC(samples, bounds_x, bounds_z)

# Interpolate random scatter onto a plotting grid
grid_x, grid_z, grid_vals = interpolate_to_grid(mc_x, mc_z, mc_status, bounds_x, bounds_z, 60)

# Separate raw data for the scatter overlay
f_x = [mc_x[i] for i in 1:samples if mc_status[i] == 1.0]
f_z = [mc_z[i] for i in 1:samples if mc_status[i] == 1.0]
inf_x = [mc_x[i] for i in 1:samples if mc_status[i] == 0.0]
inf_z = [mc_z[i] for i in 1:samples if mc_status[i] == 0.0]

# ---- Plotting ---- #
# Plot the interpolated contour boundary at the 0.5 probability threshold
boundary_plot = contour(
    grid_x, grid_z, grid_vals, 
    levels=[0.5], 
    color=:black, 
    linewidth=2,
    legend=:outerright,
    title="Controllable Region Estimation (MC + IDW)", 
    xlabel="Initial X Position [m]", 
    ylabel="Initial Z Position [m]",
    size=(800, 600)
)

scatter!(boundary_plot, f_x, f_z, label="Feasible (Sample)", color=:blue, markersize=3, alpha=0.5)
scatter!(boundary_plot, inf_x, inf_z, label="Infeasible (Sample)", color=:red, markersize=3, alpha=0.5)

display(boundary_plot)
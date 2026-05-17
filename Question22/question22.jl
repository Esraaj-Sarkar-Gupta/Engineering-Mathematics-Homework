# ===== Engineering Mathematics -- Question 21 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 2nd May, 2026

# Optimal Control with perturbations

using JuMP
using Ipopt
using Plots

g = 1.62

# ---- Define System ---- #
function system(state, u, m)
    x, z, vx, vz = state
    ux, uz = u

    xdot = vx
    zdot = vz

    vxdot = ux
    vzdot = uz - g

    dState = xdot, zdot, vxdot, vzdot
    return dState
end

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
    # No Support
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

        ux[1:N+1]
        uz[1:N+1]
    end)

    # -- Solver Constrains -- #
    for k in 1:N
        if disturbance
            α = rand() * (0.9 - 0.7) + 0.7
        else
            α = 1
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
            sqrt(ux[k]^2 + uz[k]^2 + 1e-6) * h for k in 1:N+1
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
    @constraint(model, [k=1:length(ux)], ux[k]^2 + uz[k]^2 .<= u_max^2)

    optimize!(model)

    time_array = range(t_start, t_end, length = N+1)
    return time_array, value.(x), value.(z), value.(vx), value.(vz), value.(ux), value.(uz)

end

# ---- Running it all ---- #
N = 100

# ---- System Boundary Conditions ---- #
state_zero = 5, 10, -1, -2
state_end = 0, 0, 0, 0
tspan = (0.0, 5.0)

method =:euler

t, xx, xz, vx, vz, ux, uz = optimize_trajectory(
    N,
    state_zero,
    state_end,
    tspan,
    10,
    true,
    method
)


# ---- Plotting ---- #
thrust_mag = sqrt.(ux.^2 + uz.^2)

p1 = plot(xx, xz, label=string(method), title="Trajectory (X vs Z)", xlabel="X [m]", ylabel="Z [m]", lw=2)
p2 = plot(t, xx, label=string(method), title="Trajectory (X vs t)", xlabel="t [s]", ylabel="X [m]", lw=2)
p3 = plot(t, xz, label=string(method), title="Trajectory (Z vs t)", xlabel="t [s]", ylabel="Z [m]", lw=2)

p4 = plot(t, ux, label="u_x", lw=2, title="Control Inputs vs Time", xlabel="t [s]", ylabel="Acceleration [m/s²]")
plot!(p4, t, uz, label="u_z", lw=2)
plot!(p4, t, thrust_mag, label="|u| (Total)", lw=2, linestyle=:dash, color=:black)

final_plot = plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 800))

display(final_plot)
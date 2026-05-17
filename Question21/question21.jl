# ===== Engineering Mathematics -- Question 21 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 2nd May, 2026

# Optimal Control: Euler vs Trapezoid Comparison

using JuMP
using Ipopt
using Plots

g = 1.62

println("---- Run ---- ")

function system!(state, u, m)
    # Unpack state
    X_x, X_z, v_x, v_z = state
    u_x, u_z = u

    dX_x = v_x
    dX_z = v_z

    dv_x = u_x
    dv_z = u_z - g

    dState = dX_x, dX_z, dv_x, dv_z
    return dState
end

# ---- Optimization FUnction ---- #
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
        # Declare model variables
        @variables(model, begin
            Xx[1:N+1]
            Xz[1:N+1]

            vx[1:N+1]
            vz[1:N+1]

            ux[1:N+1]
            uz[1:N+1]
        end)
        
        for k in 1:N
            @constraint(model, Xx[k+1] == Xx[k] + h * vx[k])
            @constraint(model, Xz[k+1] == Xz[k] + h * vz[k])

            @constraint(model, vx[k+1] == vx[k] + h * ux[k])
            @constraint(model, vz[k+1] == vz[k] + h * (uz[k] - g))
        end

    elseif (method == :trapezoid)
        # Declare model variables
        @variables(model, begin
            Xx[1:N+1]
            Xz[1:N+1]

            vx[1:N+1]
            vz[1:N+1]

            ux[1:N+1]
            uz[1:N+1]
        end)

        for k in 1:N
            @constraint(model, Xx[k+1] == Xx[k] + h/2 * (vx[k] + vx[k+1]))
            @constraint(model, Xz[k+1] == Xz[k] + h/2 * (vz[k] + vz[k+1]))

            @constraint(model, vx[k+1] == vx[k] + h/2 * (ux[k] + ux[k+1]))
            @constraint(model, vz[k+1] == vz[k] + h/2 * (uz[k] + uz[k+1] - 2*g))
        end
    else
        println("Optimization method was not recognized")
        return true
    end

    @objective(
        model,
        Min,
        sum(
            sqrt(ux[k]^2 + uz[k]^2 + 1e-6) * h for k in 1:N+1
        )
    )

    # Define Model Constrains
    @constraint(model, Xx[1] == state_zero[1])
    @constraint(model, Xx[N+1] == state_end[1])

    @constraint(model, Xz[1] == state_zero[2])
    @constraint(model, Xz[N+1] == state_end[2])

    @constraint(model, vx[1] == state_zero[3])
    @constraint(model, vx[N+1] == state_end[3])

    @constraint(model, vz[1] == state_zero[4])
    @constraint(model, vz[N+1] == state_end[4])

    # State Constraint
    @constraint(model, Xz .>= 0)

    # Energy Constraint: ux^2 + uz^2 <= u_max^2
    @constraint(model, [k=1:length(ux)], ux[k]^2 + uz[k]^2 <= u_max^2)

    optimize!(model)

    if (method == :euler)
        time_array = range(t_start, t_end, length = N+1)
    elseif (method == :trapezoid)
        time_array = range(t_start, t_end, length = N+1)
    end

    return time_array, value.(Xx), value.(Xz), value.(vx), value.(vz), value.(ux), value.(uz)
end

# ---- Running it all ---- #
N = 100

# ---- System Boundary Conditions ---- #
state_zero = 5, 10, -1, -2
state_end = 0, 0, 0, 0
tspan = (0.0, 5.0)

method = :trapezoid

t, xx, xz, vx, vz, ux, uz = optimize_trajectory(
    N,
    state_zero,
    state_end,
    tspan,
    10,
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
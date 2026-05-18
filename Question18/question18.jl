#
# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Dtae: 25th April, 2026

# Optimal Control

using JuMP
using Ipopt
using Plots

# ---- Define System Dynamics ---- #
function system!(state, u, m)
    a = u * m

    x, v = state

    dv = a
    dx = v

    dState = dx, dv
    return dState
end

# -- Initial Conditions -- #
t_span = (0.0, 1.0)
x_span = (0.0, 1.0)
v_span = (0.0, 0.0)

# ---- Optimization Main ---- #
# Using JuMP + Ipopt instead of fmincon

function optimize_trajectory(
    N;
    method = :euler
)
    t_start, t_end = t_span
    h = (t_end - t_start) / N

    # Initialize the Optimization Model
    model = Model(Ipopt.Optimizer)
    set_silent(model) # Shush

    # Define unknowns
    @variables(model, begin
        x[1:N+1]
        v[1:N+1]
        u[1:N+1]
    end)

    # -- Enforce Constraints -- #

    # Boundary
    @constraint(model, x[1] == x_span[1])
    @constraint(model, x[N+1] == x_span[2])
    
    @constraint(model, v[1] == v_span[1])
    @constraint(model, v[N+1] == v_span[2])

    # Physics
    if method == :euler
        for k in 1:N
            @constraint(model, x[k+1] == x[k] + h * v[k])
            @constraint(model, v[k+1] == v[k] + h * u[k])
        end

        # Constrain ghost variable
        @constraint(model, u[N+1] == u[N])

        @objective(
            model,
            Min,
            sum(u[k]^2 * h for k in 1:N)
        )

    elseif method == :trapezoid
        for k in 1:N
            @constraint(model, x[k+1] == x[k] + h/2 * (v[k] + v[k+1]))
            @constraint(model, v[k+1] == v[k] + h/2 * (u[k] + u[k+1]))
        end

        @objective(
            model,
            Min,
            sum(u[k]^2 * h for k in 1:N+1) - 
            (u[1]^2 * h + u[N+1]^2 * h)/2 # Minimize Energy
        )
    else
        error("Unrecognized integration method. Use :euler or :trapezoid")
    end

    # Solver
    optimize!(model) # Yay

    time_array = range(t_span[1], t_span[2], length = N+1)

    return time_array, value.(x), value.(v), value.(u)
end

# ---- RUNNING (and stuff) ---- #
N = 50
t_opt, x_opt, v_opt, u_opt = 
    optimize_trajectory(
        N,
        method=:euler
        )

# Plotting
p_x = plot(t_opt, x_opt, label="Position x(t)", linewidth=3, color=:blue, ylabel="x", title="Optimal Trajectory (N=$N)")
p_v = plot(t_opt, v_opt, label="Velocity v(t)", linewidth=3, color=:orange, ylabel="v")
p_u = plot(t_opt, u_opt, label="Force u(t)", linewidth=3, color=:red, ylabel="u", xlabel="Time (s)")

p_combined = plot(p_x, p_v, p_u, layout=(3,1), size=(600, 800), grid=true)
display(p_combined)




    
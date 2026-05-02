# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 25th April, 2026

# Optimal Control: Euler vs Trapezoid Comparison

using JuMP
using Ipopt
using Plots

println("Running Optimization Comparison... hold onto your matrices.")

# ---- Define System Dynamics ---- #
# (Left here for your structural comfort, even though JuMP handles the physics directly below)
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

    # Define Objective
    @objective(
        model,
        Min,
        sum(u[k]^2 * h for k in 1:N+1) # Minimize Energy
    )

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
    elseif method == :trapezoid
        for k in 1:N
            @constraint(model, x[k+1] == x[k] + h/2 * (v[k] + v[k+1]))
            @constraint(model, v[k+1] == v[k] + h/2 * (u[k] + u[k+1]))
        end
    else
        error("Unrecognized integration method. Try typing it correctly.")
    end

    # Solver
    optimize!(model) 

    time_array = range(t_span[1], t_span[2], length = N+1)

    return time_array, value.(x), value.(v), value.(u)
end

# ---- RUNNING COMPARISON ---- #
N = 20 # Large enough to be smooth, small enough to expose Euler's lies
t_euler, x_euler, v_euler, u_euler = optimize_trajectory(N, method=:euler)
t_trap, x_trap, v_trap, u_trap = optimize_trajectory(N, method=:trapezoid)

# ---- Plotting ---- #
# Position Plot
p_x = plot(
    t_euler, x_euler, 
    label="Euler", linewidth=3, color=:blue, alpha=0.5, 
    ylabel="Position (x)", title="Optimal Trajectory Comparison (N=$N)", legend=:topleft
)
plot!(
    p_x, t_trap, x_trap, 
    label="Trapezoid", linewidth=3, color=:cyan, linestyle=:dash
)

# Velocity Plot
p_v = plot(
    t_euler, v_euler, 
    label="Euler", linewidth=3, color=:orange, alpha=0.5, 
    ylabel="Velocity (v)", legend=:topleft
)
plot!(
    p_v, t_trap, v_trap, 
    label="Trapezoid", linewidth=3, color=:yellow, linestyle=:dash
)

# Control Force Plot (Where the magic/cheating happens)
p_u = plot(
    t_euler, u_euler, 
    label="Euler", linewidth=3, color=:red, alpha=0.5, 
    ylabel="Force (u)", xlabel="Time (s)", legend=:topright
)
plot!(
    p_u, t_trap, u_trap, 
    label="Trapezoid", linewidth=3, color=:pink, linestyle=:dash
)

# Stack 'em up
p_combined = plot(p_x, p_v, p_u, layout=(3,1), size=(800, 1000), grid=true)
display(p_combined)
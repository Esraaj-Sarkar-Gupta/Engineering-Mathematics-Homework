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

    # Define Objective
    @objective(
        model,
        Min,
        sum(u[k]^2 * h for k in 1:N) # Minimize Energy
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
        error("Unrecognized integration method. Use :euler or :trapezoid")
    end

    # Solver
    optimize!(model) # Yay

    time_array = range(t_span[1], t_span[2], length = N+1)

    return time_array, value.(x), value.(v), value.(u)
end

# ---- RUNNING (and stuff) ---- #
N = 5
t_opt, x_opt, v_opt, u_opt = optimize_trajectory(N, method=:trapezoid)

# Plotting
p_x = plot(t_opt, x_opt, label="Position x(t)", linewidth=3, color=:blue, ylabel="x", title="Optimal Trajectory (N=$N)")
p_v = plot(t_opt, v_opt, label="Velocity v(t)", linewidth=3, color=:orange, ylabel="v")
p_u = plot(t_opt, u_opt, label="Force u(t)", linewidth=3, color=:red, ylabel="u", xlabel="Time (s)")

p_combined = plot(p_x, p_v, p_u, layout=(3,1), size=(600, 800), grid=true)
display(p_combined)

# ---- Animation: The Evolution of Euler's Incompetence ---- #

println("Fetching the popcorn... and a fire extinguisher for your CPU. Generating animation.")

lims_x = (-0.1, 1.1)
lims_v = (-0.2, 1.8)
lims_u = (-10, 10)

# Create the animation object
anim = @animate for N in 2:2:100
    t_euler, x_euler, v_euler, u_euler = optimize_trajectory(N, method=:euler)
    t_trap, x_trap, v_trap, u_trap = optimize_trajectory(N, method=:trapezoid)
    
    # Position
    p_x = plot(
        t_euler, x_euler, 
        label="Euler", lw=4, color=:blue, alpha=0.4, 
        ylabel="Position (x)", ylims=lims_x, legend=:topleft
    )
    plot!(
        p_x, t_trap, x_trap, 
        label="Trapezoid", lw=2, color=:black, linestyle=:dash
    )
    
    # Velocity
    p_v = plot(
        t_euler, v_euler, 
        label="Euler", lw=4, color=:orange, alpha=0.4, 
        ylabel="Velocity (v)", ylims=lims_v, legend=:topleft
    )
    plot!(
        p_v, t_trap, v_trap, 
        label="Trapezoid", lw=2, color=:black, linestyle=:dash
    )
    
    # Control Force
    p_u = plot(
        t_euler, u_euler, 
        label="Euler", lw=4, color=:red, alpha=0.4, 
        ylabel="Force (u)", xlabel="Time (s)", ylims=lims_u, legend=:topright
    )
    plot!(
        p_u, t_trap, u_trap, 
        label="Trapezoid", lw=2, color=:black, linestyle=:dash
    )
    
    plot(
        p_x, p_v, p_u, 
        layout=(3,1), 
        size=(800, 1000), 
        grid=true,
        plot_title="Optimal Control Comparison: N = $N"
    )
end

# Save
gif_filename = "euler_vs_trapezoid_evolution.gif"
gif(anim, gif_filename, fps = 10)

println("Blockbuster hit saved as $gif_filename. Check your local directory.")
# ===== Engineering Mathematics -- Homework 5 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 5th April, 2026

# This script generates the periodic solutions simple second order systems.

using Plots
using DifferentialEquations

# ---- General Second Order Systems ---- #

function second_order_system!(du, u, p, t)
    x, y, vx, vy = u
    k, m = p

    K₁ = -k/m
    K₂ = 0.2 * K₁

    # Physics
    # r̈ = - k/m

    # Define system dynamics matrix
    A = [
        0 0 1 0;
        0 0 0 1;
        K₁ 0 0 0;
        0 K₁ 0 0
    ]

    du .= A * u
end

tspan = (0, 50)
u₀ = [1, 1, 1, 1,]

p = (
    10, # k
    5,  # m
)

# Fifth Order Solver
dynamics_problem = ODEProblem(second_order_system!,u₀, tspan, p)
solution = solve(dynamics_problem, Tsit5(), saveat=0.01)

# Extract solutions
x_sol = solution[1, :]
y_sol = solution[2, :]
times = solution.t

# Plot
p1 = plot(
    x_sol,
    y_sol,
    xlabel = "x",
    ylabel = "y",
    title = "Physical Space X - Y",
    color = :blue,
    grid = true,
    legend = false
)

p2 = plot(
    x_sol,
    times,
    xlabel = "Time",
    ylabel = "r",
    title = "r vs Time",
    color = :green,
    grid = true,
    label = "X"
)

plot!(
    p2,
    y_sol,
    times,
    color = :red,
    label = "Y"
)

p12 = plot(p1, p2, layout=(1,2))
png(p12, "graph.png")

# The lines used to animate the system above have been commented
# out to save compute at runtime.
"""
# Animate solution
xlims = (minimum(x_sol) - 0.5, maximum(x_sol) + 0.5)
ylims = (minimum(y_sol) - 0.5, maximum(y_sol) + 0.5)

animation = @animate for i in 1:10:length(times)
    p1 = plot(
        x_sol[1:i],
        y_sol[1:i],
        xlabel = "x",
        ylabel = "y",
        title = "Physical Space X - Y",
        color = :blue,
        grid = true,
        legend = false,
        xlims = xlims,
        ylims = ylims
    )

    scatter!(
        p1,
        [x_sol[i]],
        [y_sol[i]],
        color = :red
    )
    plot(p1, size=(800, 400))
end
"""
#gif(animation, "general_coupled_oscillator.gif", fps = 30)
#println("Animation saved!")

# ---- Inverse Square Central Gravity System ---- #
"""
function central_gravity!(du, u, p, t)
    x, y, vx, vy = u

    # For simplicity, here we consider G/M = 10

    r_sq = x^2 + y^2
    r_cu = 10 * r_sq^(1.5)

    du .= [vx, vy, -x/r_cu, -y/r_cu]
end

gu₀ = [0.1,0.1,0.5,0]
tspang = (0, 10)

ISCG_problem = ODEProblem(central_gravity!, gu₀, tspang, 0)
solution_g = solve(ISCG_problem, Tsit5(), saveat=0.001)

p3 = plot(
    solution_g[1,:],
    solution_g[2,:],
    color = :blue,
    legend = true,
    grid = true,
    xlabel = "X",
    ylabel = "Y",
    title = "Inverse Squared Central Gravity Problem",
    label = "Trajectory"
)

# Plot the attractor
scatter!(
    p3,
    [0],
    [0],
    color=:red,
    label = "Attractor"
)

# Plot the starting point
scatter!(
    p3,    
    [gu₀[1]],
    [gu₀[2]],
    color=:green,
    label = "Initial Position"
)

plot(p3)

"""

# ---- Generalising the Central Gravity System ---- #
"""
function second_order_nonlinear_system!(du, u, p, t)
    x, y, vx, vy = u
    kx, ky, n = p
    # Define System dynamics
    # Cannot use matrices since this system is non_linear

    dx = vx
    dy = vy

    dvx = - kx * x^(n-1) * x
    dvy = - ky * y^(n-1) * y

    du .= [dx, dy, dvx, dvy]
end

# Define problem parameters
tspan = (0, 10)
u₀ = [1, 0, 0, 1]

p = (
    2,  # kx
    2,  # ky
    3,  # n
)

non_linear_problem = ODEProblem(second_order_nonlinear_system!, u₀, tspan, p)
solution = solve(non_linear_problem, Tsit5(), saveat=0.01)

# Plot
p3_5 = plot(
    solution[1,:],
    solution[2,:],
    color=:magenta,
    grid=true,
    legend = false,
    title = "General Non-Linear Second Order System",
    xlabel = "X",
    ylabel = "Y",
    aspect_ratio=:equal,
)

png(p3_5, "nonLinear_secondOrder.png")

xlims = (minimum(solution[1,:]) - 0.5, maximum(solution[1,:]) + 0.5)
ylims = (minimum(solution[2,:]) - 0.5, maximum(solution[2,:]) + 0.5)

# Animate solutions
animation = @animate for i in 1:10:length(solution.t)
    p4 = plot(
        solution[1,1:i],
        solution[2,1:i],
        color=:cyan,
        grid = true,
        legend = false,
        title = "General Non-Linear Second Order System",
        xlabel = "X",
        ylabel = "Y",
        xlims = xlims,
        ylims = ylims,
        aspect_ratio=:equal,
    )

    plot!(
        p4,
        [0],
        [0],
        color=:red,
    )
end

gif(animation, "nonLinear_secondOrder.gif", fps = 30)
print("Saved Animation")
"""
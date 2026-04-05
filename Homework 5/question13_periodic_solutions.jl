# ===== Engineering Mathematics -- Homework 5 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 5th April, 2026

# This script generates the periodic solutions simple second order systems.

using Plots
using DifferentialEquations

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
        0 K₂ 0 0
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

plot(p1, p2, layout=(1,2))




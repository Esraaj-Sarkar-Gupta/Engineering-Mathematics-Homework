# ===== Engineering Mathematics -- Homework 5 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 5th April, 2026

# Inverse Square Central Gravity System 

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
solution_g = solve(ISCG_problem, Tsit5(), saveat=0.001, reltol=1e-8, abstol=1e-8)

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

display(p3)


# ---- Generalising the Central Gravity System ---- #

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
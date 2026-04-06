# ===== Engineering Mathematics -- Homework 6 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 6th April, 2026

# This script is used to generate an Unstable Osicalltor

using Plots
using DifferentialEquations

# ---- Lorenz Attractor System ---- #
function lorenz_attractor!(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p

    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

# Setup and Solve
lu₀ = [1.0, 1.0, 1.0]
p_l = (10.0, 28.0, 8/3)

tspan_pre = (0.0, 10.0)
prob_pre = ODEProblem(lorenz_attractor!, lu₀, tspan_pre, p_l)
solution_pre = solve(prob_pre, Tsit5(), saveat=0.01)

lu₁ = copy(solution_pre.u[end])
lu₁[1] += 50.0 
lu₁[2] += -5.0 
lu₁[3] += 5.0

tspan_post = (10.0, 30.0)
prob_post = ODEProblem(lorenz_attractor!, lu₁, tspan_post, p_l)
solution_post = solve(prob_post, Tsit5(), saveat=0.01)

times = vcat(solution_pre.t, solution_post.t)
x_sol = vcat(solution_pre[1,:], solution_post[1,:])
y_sol = vcat(solution_pre[2,:], solution_post[2,:])
z_sol = vcat(solution_pre[3,:], solution_post[3,:])

prob_full = ODEProblem(lorenz_attractor!, lu₀, (0.0, 30.0), p_l)
solution_full = solve(prob_full, Tsit5(), saveat=0.01)
x_unp = solution_full[1,:]
y_unp = solution_full[2,:]
z_unp = solution_full[3,:]

px_perturb = solution_pre[1,end]
py_perturb = solution_pre[2,end]
pz_perturb = solution_pre[3,end]

px_post = lu₁[1]
py_post = lu₁[2]
pz_post = lu₁[3]

xlims = (min(minimum(x_sol), minimum(x_unp)) - 5.0, max(maximum(x_sol), maximum(x_unp)) + 5.0)
ylims = (min(minimum(y_sol), minimum(y_unp)) - 5.0, max(maximum(y_sol), maximum(y_unp)) + 5.0)
zlims = (min(minimum(z_sol), minimum(z_unp)) - 5.0, max(maximum(z_sol), maximum(z_unp)) + 5.0)

idx_10 = searchsortedfirst(times, 10.0)

animation = @animate for i in 1:10:length(times)
    p1 = plot(
        x_unp[1:i],
        y_unp[1:i],
        z_unp[1:i],
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",
        title = "Perturbed Lorenz Attractor",
        color = :cyan,
        grid = true,
        legend = false,
        xlims = xlims,
        ylims = ylims,
        zlims = zlims
    )

    if i > idx_10
        plot!(
            p1,
            x_sol[idx_10:i],
            y_sol[idx_10:i],
            z_sol[idx_10:i],
            color = :blue
        )

        plot!(
            p1,
            [px_perturb, px_post],
            [py_perturb, py_post],
            [pz_perturb, pz_post],
            arrow = true,
            color = :orange,
            linewidth = 2
        )

        scatter!(
            p1,
            [px_perturb],
            [py_perturb],
            [pz_perturb],
            color = :orange,
            markersize = 6,
            marker = :star
        )
    end

    scatter!(
        p1,
        [x_unp[i]],
        [y_unp[i]],
        [z_unp[i]],
        color = :cyan,
        markersize = 3
    )

    scatter!(
        p1,
        [x_sol[i]],
        [y_sol[i]],
        [z_sol[i]],
        color = :red,
        markersize = 4
    )
end

gif(animation, "perturbed_lorenz.gif", fps = 30)
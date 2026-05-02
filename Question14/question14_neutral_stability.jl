# ===== Engineering Mathematics -- Homework 6 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 6th April, 2026

# This script is used to generate a Neutrally Stable Osicalltor

using DifferentialEquations
using Plots

# Simplest possible system for this -- Harmonic Osicalltor
# ---- Harmonic Oscillator System ---- #

function harmonic_oscillator!(du, u, p, t)
    x, y, vx, vy = u
    kx, ky = p

    du[1] = vx
    du[2] = vy
    du[3] = -kx * x
    du[4] = -ky * y
end

hu₀ = [1, 1, 0, 0]
tspan_h = (0, 20)

p_h = (
    1, 
    2, 
)

harmonic_problem = ODEProblem(harmonic_oscillator!, hu₀, tspan_h, p_h)
solution_h = solve(harmonic_problem, Tsit5(), saveat=0.01)

p1 = plot(
    solution_h[1,:],
    solution_h[2,:],
    color = :blue,
    legend = false,
    grid = true,
    xlabel = "X",
    ylabel = "Y",
    title = "Harmonic Oscillator Problem"
)

plot(p1)

# ===== Engineering Mathematics -- Homework 6 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 6th April, 2026

# This script is used to generate a Neutrally Stable Osicalltor

using DifferentialEquations
using Plots

# Simplest possible system for this -- Harmonic Osicalltor
# ---- Harmonic Oscillator System ---- #

function harmonic_oscillator!(du, u, p, t)
    x, y, vx, vy = u
    k = p

    du[1] = vx
    du[2] = vy
    du[3] = -k * x
    du[4] = -k * y
end

hu₀ = [1, 1, 0, 0]
tspan_h = (0, 20)

p_h = 2.0

harmonic_problem = ODEProblem(harmonic_oscillator!, hu₀, tspan_h, p_h)
solution_h = solve(harmonic_problem, Tsit5(), saveat=0.01)

p1 = plot(
    solution_h[1,:],
    solution_h[2,:],
    color = :blue,
    legend = false,
    grid = true,
    xlabel = "X",
    ylabel = "Y",
    title = "Harmonic Oscillator Problem"
)

plot(p1)

# ---- Animation (with perturbance) ---- #
hu₀ = [1.0, 1.0, 0.0, 0.0]
p_h = 2.0

tspan_pre = (0.0, 10.0)
prob_pre = ODEProblem(harmonic_oscillator!, hu₀, tspan_pre, p_h)
solution_pre = solve(prob_pre, Tsit5(), saveat=0.01)

hu₁ = copy(solution_pre.u[end])
hu₁[3] += 2.0 
hu₁[4] += -1.5 

tspan_post = (10.0, 30.0)
prob_post = ODEProblem(harmonic_oscillator!, hu₁, tspan_post, p_h)
solution_post = solve(prob_post, Tsit5(), saveat=0.01)

times = vcat(solution_pre.t, solution_post.t)
x_sol = vcat(solution_pre[1,:], solution_post[1,:])
y_sol = vcat(solution_pre[2,:], solution_post[2,:])

px_perturb = solution_pre[1,end]
py_perturb = solution_pre[2,end]

xlims = (minimum(x_sol) - 0.5, maximum(x_sol) + 0.5)
ylims = (minimum(y_sol) - 0.5, maximum(y_sol) + 0.5)

idx_10 = searchsortedfirst(times, 10.0)

animation = @animate for i in 1:10:length(times)
    p1 = plot(
        x_sol[1:min(i, idx_10)],
        y_sol[1:min(i, idx_10)],
        xlabel = "X",
        ylabel = "Y",
        title = "Perturbed Harmonic Oscillator",
        color = :cyan,
        grid = true,
        legend = false,
        xlims = xlims,
        ylims = ylims
    )

    if i > idx_10
        plot!(
            p1,
            x_sol[idx_10:i],
            y_sol[idx_10:i],
            color = :blue
        )

        plot!(
            p1,
            [px_perturb, px_perturb + 1.0],
            [py_perturb, py_perturb - 0.75],
            arrow = true,
            color = :orange,
            linewidth = 2
        )

        scatter!(
            p1,
            [px_perturb],
            [py_perturb],
            color = :orange,
            markersize = 8,
            marker = :star
        )
    end

    scatter!(
        p1,
        [x_sol[i]],
        [y_sol[i]],
        color = :red
    )
end

gif(animation, "perturbed_oscillator.gif", fps = 30)

# ---- Save last frame of animation ---- #

p2 = plot(
        x_sol[1:idx_10],
        y_sol[1:idx_10],
        xlabel = "X",
        ylabel = "Y",
        title = "Perturbed Harmonic Oscillator",
        color = :cyan,
        grid = true,
        legend = false,
        xlims = xlims,
        ylims = ylims
    )

    if true
        plot!(
            p2,
            x_sol[idx_10:length(x_sol)],
            y_sol[idx_10:length(y_sol)],
            color = :blue
        )

        plot!(
            p2,
            [px_perturb, px_perturb + 1.0],
            [py_perturb, py_perturb - 0.75],
            arrow = true,
            color = :orange,
            linewidth = 2
        )

        scatter!(
            p2,
            [px_perturb],
            [py_perturb],
            color = :orange,
            markersize = 8,
            marker = :star
        )
    end

plot(p2)
png(p2, "perturbed_harmonic_oscillator.png")
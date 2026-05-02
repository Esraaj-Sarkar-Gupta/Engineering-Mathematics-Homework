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

p_final = plot(
    x_unp,
    y_unp,
    z_unp,
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

plot!(
    p_final,
    x_sol[idx_10:end],
    y_sol[idx_10:end],
    z_sol[idx_10:end],
    color = :blue
)

plot!(
    p_final,
    [px_perturb, px_post],
    [py_perturb, py_post],
    [pz_perturb, pz_post],
    arrow = true,
    color = :orange,
    linewidth = 2
)

scatter!(
    p_final,
    [px_perturb],
    [py_perturb],
    [pz_perturb],
    color = :orange,
    markersize = 6,
    marker = :star
)

scatter!(
    p_final,
    [x_unp[end]],
    [y_unp[end]],
    [z_unp[end]],
    color = :cyan,
    markersize = 3
)

scatter!(
    p_final,
    [x_sol[end]],
    [y_sol[end]],
    [z_sol[end]],
    color = :red,
    markersize = 4
)

savefig(p_final, "perturbed_lorenz_final.png")

gif(animation, "perturbed_lorenz.gif", fps = 30)


# ---- Unstable Circular Repeller System ---- #

function unstable_circle!(du, u, p, t)
    x, y = u
    R = p[1]

    r_sq = x^2 + y^2

    du[1] = y + x * (r_sq - R^2)
    du[2] = -x + y * (r_sq - R^2)
end

u₀_perf = [1.0, 0.0]
tspan_perf = (0.0, 15.0)
p_c = (1.0,)

prob_perf = ODEProblem(unstable_circle!, u₀_perf, tspan_perf, p_c)
sol_perf = solve(prob_perf, Tsit5(), saveat=0.01)

u₀_in = [0.99, 0.0]
tspan_in = (0.0, 15.0)
prob_in = ODEProblem(unstable_circle!, u₀_in, tspan_in, p_c)
sol_in = solve(prob_in, Tsit5(), saveat=0.01)

u₀_out = [1.01, 0.0]
tspan_out = (0.0, 3.5)
prob_out = ODEProblem(unstable_circle!, u₀_out, tspan_out, p_c)
sol_out = solve(prob_out, Tsit5(), saveat=0.01)

p_c_plot = plot(
    sol_perf[1,:],
    sol_perf[2,:],
    color = :black,
    linewidth = 2,
    legend = true,
    grid = true,
    xlabel = "X",
    ylabel = "Y",
    title = "Unstable Circular Limit Cycle",
    label = "Unstable Orbit (R=1)"
)

plot!(
    p_c_plot,
    sol_in[1,:],
    sol_in[2,:],
    color = :blue,
    label = "Starts Inside (Collapses)"
)

plot!(
    p_c_plot,
    sol_out[1,:],
    sol_out[2,:],
    color = :red,
    label = "Starts Outside (Escapes)"
)

plot(p_c_plot)

# ---- Reverse van der Pol Osicalltor ---- #


using Plots
using DifferentialEquations

function vdp_normal!(du, u, p, t)
    x, vx = u
    μ = p[1]
    
    du[1] = vx
    du[2] = μ * (1.0 - x^2) * vx - x
end

function vdp_reversed!(du, u, p, t)
    x, vx = u
    μ = p[1]
    
    du[1] = vx
    du[2] = -μ * (1.0 - x^2) * vx - x
end

p_vdp = (1.5,)

u0_ref = [0.1, 0.0]
tspan_ref = (0.0, 50.0)
prob_ref = ODEProblem(vdp_normal!, u0_ref, tspan_ref, p_vdp)
sol_ref = solve(prob_ref, Tsit5(), saveat=0.01)

idx_tail = searchsortedfirst(sol_ref.t, 40.0)

u0_in = [1.0, 1.0]
tspan_in = (0.0, 20.0)
prob_in = ODEProblem(vdp_reversed!, u0_in, tspan_in, p_vdp)
sol_in = solve(prob_in, Tsit5(), saveat=0.01)

u0_out = [2.5, 2.5]
tspan_out = (0.0, 5.0)
prob_out = ODEProblem(vdp_reversed!, u0_out, tspan_out, p_vdp)
sol_out = solve(prob_out, Tsit5(), saveat=0.01)

p_vdp_plot = plot(
    sol_ref[1,idx_tail:end],
    sol_ref[2,idx_tail:end],
    color = :black,
    linewidth = 2,
    linestyle = :dash,
    legend = true,
    grid = true,
    xlabel = "Position (X)",
    ylabel = "Velocity (V)",
    title = "Reversed Van der Pol (Unstable Limit Cycle)",
    label = "Unstable Orbit (Repeller)"
)

plot!(
    p_vdp_plot,
    sol_in[1,:],
    sol_in[2,:],
    color = :blue,
    label = "Starts Inside (Collapses to Origin)"
)

plot!(
    p_vdp_plot,
    sol_out[1,:],
    sol_out[2,:],
    color = :red,
    label = "Starts Outside (Escapes)"
)

plot(p_vdp_plot)
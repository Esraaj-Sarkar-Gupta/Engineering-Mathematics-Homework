using Plots
using DifferentialEquations

"""
function vdp_stable!(du, u, p, t)
    x, vx = u
    μ = p[1]
    
    du[1] = vx
    du[2] = μ * (1.0 - x^2) * vx - x
end

p_vdp = (1.5,)
tspan = (0.0, 300.0)

u0_outside = [3.0, 3.0]
prob_out = ODEProblem(vdp_stable!, u0_outside, tspan, p_vdp)
sol_out = solve(prob_out, Tsit5(), saveat=0.01)

u0_inside = [0.1, 0.1]
prob_in = ODEProblem(vdp_stable!, u0_inside, tspan, p_vdp)
sol_in = solve(prob_in, Tsit5(), saveat=0.01)

p_stable = plot(
    sol_out[1, end-2000:end],
    sol_out[2, end-2000:end],
    color = :black,
    linewidth = 4,
    alpha = 0.3, # Lower alpha makes it a background track
    label = "Stable Limit Cycle",
    grid = true,
    xlabel = "Position (X)",
    ylabel = "Velocity (V)",
    title = "Asymptotically Stable Limit Cycle",
    xlims = (-3, 3.5),
    ylims = (-3.5, 3.5)
)

plot!(
    p_stable,
    sol_out[1,:],
    sol_out[2,:],
    color = :red,
    alpha = 0.8, # Increased opacity to see it over the black line
    linewidth = 1.5,
    label = "Outside Start (Spirals In)"
)

plot!(
    p_stable,
    sol_in[1,:],
    sol_in[2,:],
    color = :blue,
    alpha = 0.8, # Increased opacity
    linewidth = 1.5,
    label = "Inside Start (Spirals Out)"
)

plot(p_stable)
"""

# Animation

function vdp_stable!(du, u, p, t)
    x, vx = u
    μ = p[1]
    
    du[1] = vx
    du[2] = μ * (1.0 - x^2) * vx - x
end

p_vdp = (1.5,)
tspan = (0.0, 100.0)

u0_outside = [3.0, 3.0]
prob_out = ODEProblem(vdp_stable!, u0_outside, tspan, p_vdp)
sol_out = solve(prob_out, Tsit5(), saveat=0.1)

u0_inside = [0.1, 0.1]
prob_in = ODEProblem(vdp_stable!, u0_inside, tspan, p_vdp)
sol_in = solve(prob_in, Tsit5(), saveat=0.1)

x_lc = sol_out[1, end-200:end]
v_lc = sol_out[2, end-200:end]

animation = @animate for i in 1:length(sol_out.t)
    p = plot(
        x_lc, 
        v_lc, 
        color = :black, 
        alpha = 0.2, 
        linewidth = 3, 
        label = "Limit Cycle",
        xlabel = "Position (X)",
        ylabel = "Velocity (V)",
        title = "Van der Pol Convergence (t = $(round(sol_out.t[i], digits=1)))",
        xlims = (-3, 3.5),
        ylims = (-4, 4),
        grid = true,
        legend = :topright
    )

    plot!(
        p, 
        sol_out[1, 1:i], 
        sol_out[2, 1:i], 
        color = :red, 
        alpha = 0.8, 
        label = "Spiral In"
    )

    plot!(
        p, 
        sol_in[1, 1:i], 
        sol_in[2, 1:i], 
        color = :blue, 
        alpha = 0.8, 
        label = "Spiral Out"
    )

    scatter!(
        p, 
        [sol_out[1, i], sol_in[1, i]], 
        [sol_out[2, i], sol_in[2, i]], 
        color = [:red, :blue], 
        markersize = 4, 
        label = false
    )
end

gif(animation, "vdp_convergence.gif", fps = 30)
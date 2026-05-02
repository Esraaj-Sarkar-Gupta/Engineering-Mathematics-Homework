# ===== Engineering Mathematics -- Homework 5, Question 13 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 12th April, 2026

# Finding periodic solutions using root-finding

using DifferentialEquations
using LinearAlgebra
using Plots

using NLsolve

# ---- Fundamental Natural Constants ---- #

g = 9.810

# ---- Fundamental System Parameters ---- #

struct fparameters
    m  :: Float64
    rₐ :: Vector{Float64}
    rᵦ :: Vector{Float64}

    Lₐ :: Float64 
    Lᵦ :: Float64 
    
    kₐ :: Float64
    kᵦ :: Float64

    cₐ :: Float64
    cᵦ :: Float64

    c_d :: Float64
end

# ---- Derived System Parameters ---- #

struct dparameters
    lₐ :: Float64
    lᵦ :: Float64

    a :: Vector{Float64}
    b :: Vector{Float64}

end

function get_dparameters(sys::fparameters, r_c::Vector{Float64})
    a_vec = sys.rₐ - r_c
    b_vec = sys.rᵦ - r_c
    
    # Current scalar lengths
    lₐ_val = norm(a_vec)
    lᵦ_val = norm(b_vec)
    
    return dparameters(lₐ_val, lᵦ_val, a_vec, b_vec)
end

# ---- Define System Evolution ---- #

function spring_damper_mass_system!(du, u, p, t)
    # Unpack state
    r_c = u[1:2]
    v_c = u[3:4]

    # Unpack Parameters
    m = p.m 
    r_A, r_B = p.rₐ, p.rᵦ

    k_A, k_B = p.kₐ, p.kᵦ
    L_A, L_B = p.Lₐ, p.Lᵦ

    c_A, c_B, c_d = p.cₐ, p.cᵦ, p.c_d

    # -- Geometry -- #
    vec_to_A = r_A - r_c
    vec_to_B = r_B - r_c

    dist_A = norm(vec_to_A)
    dist_B = norm(vec_to_B)

    # Unit Vectors
    hat_A = vec_to_A / dist_A
    hat_B = vec_to_B / dist_B

    # -- Force Calculations -- #

    # - Springs - #
    delta_A = dist_A - L_A
    delta_B = dist_B - L_B
    
    F_spring_A = k_A * delta_A * hat_A
    F_spring_B = k_B * delta_B * hat_B

    # - Dampers - #
    vel_proj_A = dot(v_c, -hat_A)
    vel_proj_B = dot(v_c, -hat_B)

    F_damp_A = -c_A * dot(v_c, hat_A) * hat_A 
    F_damp_B = -c_B * dot(v_c, hat_B) * hat_B
    
    # - Drag - #
    F_drag = -c_d * v_c

    # - Gravity - #
    F_gravity = [0.0, -m * g]

    # -- Newton's Second Law -- #

    F_total = F_spring_A + F_spring_B + F_damp_A + F_damp_B + F_drag + F_gravity
    accel = F_total / m

    # -- Update State -- #
    du[1:2] = v_c    # dr/dt = v
    du[3:4] = accel  # dv/dt = F/m

end

# ---- Solving the Unremarkable System ---- #

params = fparameters(
    5.0,          # Mass (kg)
    [-3.0, 3.0],  # Anchor A
    [ 3.0,-3.0],  # Anchor B
    5.0,          # Rest Length A
    5.0,          # Rest Length B
    1.0,          # k_A
    5.0,          # k_B
    0.0,          # c_A
    0.0,          # c_B
    0.0           # c_d (Air Drag)
)

# Time Span
tspan = (0.0, 50.0)

# -- Convinience -- #

show_path = true

rA = params.rₐ
rB = params.rᵦ

LA = params.Lₐ * rA / norm(rA)

# ---- Looking for a Periodic Solution ---- #
# Using the Shooting Function Method
x_fixed = 0.0 # Poincare Manifold

function shooting_residual!(F, vars)
    # Unpack guesses
    y_guess, vx_guess, vy_guess, T_guess = vars

    # Force x = 0.0 (Poincare Manifold)
    u_start = [x_fixed, y_guess, vx_guess, vy_guess]

    # Tight tolerance solving
    prob_shoot = ODEProblem(spring_damper_mass_system!, u_start, (0.0, T_guess), params)
    sol_shoot = solve(prob_shoot, Tsit5(), reltol=1e-11, abstol=1e-11)

    u_end = sol_shoot.u[end]

    # Check residuals
    F[1] = u_end[1] - x_fixed
    F[2] = u_end[2] - y_guess
    F[3] = u_end[3] - vx_guess
    F[4] = u_end[4] - vy_guess
end

# ---- Using NLsolve + MonteCarlo ---- #
max_attempts = 100
found_solution = false  # flag
best_root = zeros(4)

# -- Helpers -- #
struct SearchBounds
    ymin::Float64
    ymax::Float64
    
    vxmin::Float64
    vxmax::Float64
    
    vymin::Float64
    vymax::Float64
    
    Tmin::Float64
    Tmax::Float64
end

function generate_random_guess(bounds::SearchBounds)
    ymin, ymax   = bounds.ymin, bounds.ymax
    vxmin, vxmax = bounds.vxmin, bounds.vxmax
    vymin, vymax = bounds.vymin, bounds.vymax
    Tmin, Tmax   = bounds.Tmin, bounds.Tmax

    y_rand  = ymin  + (ymax - ymin)   * rand()
    vx_rand = vxmin + (vxmax - vxmin) * rand()
    vy_rand = vymin + (vymax - vymin) * rand()
    T_rand  = Tmin  + (Tmax - Tmin)   * rand()

    # 3. Return the formatted guess vector
    return [y_rand, vx_rand, vy_rand, T_rand]
end

search_space = SearchBounds(
    3.0, 7.0,   # Y position bounds
    5.0, 7.0,   # Vx bounds
    4.0, 15.0,  # Vy bounds
    15.0, 20.0  # Period bounds
)

println("Starting Random Multi-Start Shooting Method...")

for attempt in 1:max_attempts
    
    # Generate the bounded random guess using our new function
    initial_guess_pos = generate_random_guess(search_space)
    
    print("\rAttempt $attempt/$max_attempts... ") 

    result_pos = try
        nlsolve(
            shooting_residual!,
            initial_guess_pos,
            autodiff=:forward,
            show_trace=false, 
            iterations=50     
        )
    catch e
        nothing 
    end

    if result_pos !== nothing && converged(result_pos) && result_pos.zero[4] > 0
        global best_root = result_pos.zero
        global found_solution = true
        println("\n\nSuccess on attempt $(attempt)!")
        break
    end
end

# ---- 3D -> 2D Projected Helical Spring Generator ---- #
function get_spring_coords(p1, p2; R=0.35, coils=12, points=200, tilt=0.6)
    v = p2 - p1
    L = norm(v)
    
    # Unit vectors (axial and normal)
    u = v / (L + 1e-6)
    n = [-u[2], u[1]]
    
    xs = zeros(points)
    ys = zeros(points)
    
    t_vals = range(0, 1, length=points)
    for (i, t) in enumerate(t_vals)
        # Angle parameter for the helix
        theta = coils * 2π * t
        
        # Taper envelope to make ends connect cleanly to the points
        env = sin(π * t) 
        
        # Calculate offset in 3D and project to 2D
        # tilt controls how much the out-of-plane loop overlaps in 2D
        u_dist = L * t + (R * tilt * env) * sin(theta)
        n_dist = (R * env) * cos(theta)
        
        pos = p1 .+ u_dist .* u .+ n_dist .* n
        xs[i] = pos[1]
        ys[i] = pos[2]
    end
    
    return xs, ys
end

# -- Representation -- #
if found_solution
    println("  x₀  = $(x_fixed) m (Fixed Manifold)")
    println("  y₀  = $(round(best_root[1], digits=4)) m")
    println("  vx₀ = $(round(best_root[2], digits=4)) m/s")
    println("  vy₀ = $(round(best_root[3], digits=4)) m/s")
    println("  T   = $(round(best_root[4], digits=4)) s")

    u_perfect = [x_fixed, best_root[1], best_root[2], best_root[3]]
    T_perfect = best_root[4]
    
    prob_perfect = ODEProblem(spring_damper_mass_system!, u_perfect, (0.0, T_perfect), params)
    sol_perfect = solve(prob_perfect, Tsit5(), saveat=0.01)

    x_perf = sol_perfect[1, :]
    y_perf = sol_perfect[2, :]

    pad_p = 2.0
    xmin_p = min(minimum(x_perf), params.rₐ[1], params.rᵦ[1]) - pad_p
    xmax_p = max(maximum(x_perf), params.rₐ[1], params.rᵦ[1]) + pad_p
    ymin_p = min(minimum(y_perf), params.rₐ[2], params.rᵦ[2]) - pad_p
    ymax_p = max(maximum(y_perf), params.rₐ[2], params.rᵦ[2]) + pad_p

    p_periodic = plot(
        x_perf, y_perf, 
        title="Fundamental Periodic Orbit (T = $(round(T_perfect, digits=2))s)", 
        xlabel="X (m)", ylabel="Y (m)",
        aspect_ratio=:equal,
        grid=true, legend=false,
        color=:purple, linewidth=2,
        xlim=(xmin_p, xmax_p), ylim=(ymin_p, ymax_p)
    )
    
    scatter!(p_periodic, [params.rₐ[1], params.rᵦ[1]], [params.rₐ[2], params.rᵦ[2]], color=:black, markersize=6)
    display(p_periodic)
    png(p_periodic, "periodic_orbit_static.png")

    println("\nGenerating 3D-projected animation for the periodic orbit...")
    
    anim_periodic = @animate for i in 1:4:length(sol_perfect.t)
        xc, yc = x_perf[i], y_perf[i]

        p = plot(
            xlim=(xmin_p, xmax_p), ylim=(ymin_p, ymax_p),
            aspect_ratio=:equal,
            xlabel="X (m)", ylabel="Y (m)",
            title="Periodic Orbit    t = $(round(sol_perfect.t[i], digits=2)) s",
            legend=false, grid=true
        )

        plot!(p, x_perf, y_perf, lw=2, alpha=0.2, color=:purple)
        
        plot!(p, x_perf[1:i], y_perf[1:i], lw=2, color=:purple)

        # Draw 3D Projected Springs
        sx_A, sy_A = get_spring_coords(params.rₐ, [xc, yc]; R=0.4, coils=12, points=300, tilt=0.5)
        sx_B, sy_B = get_spring_coords(params.rᵦ, [xc, yc]; R=0.4, coils=12, points=300, tilt=0.5)
        
        plot!(p, sx_A, sy_A, lw=1.5, color=:gray)
        plot!(p, sx_B, sy_B, lw=1.5, color=:gray)

        # Static Anchors
        scatter!(p,
            [params.rₐ[1], params.rᵦ[1]],
            [params.rₐ[2], params.rᵦ[2]],
            markersize=8, color=:black, markerstrokewidth=0
        )

        # Dynamic Mass point
        scatter!(p, [xc], [yc], markersize=10, color=:red, markerstrokewidth=0)
    end

    gif(anim_periodic, "periodic_orbit_3D_springs.gif", fps=30)
    println("Saved periodic_orbit_3D_springs.gif!")
else
    println("\nFailed to find a periodic orbit after $max_attempts random attempts.")
end
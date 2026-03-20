# ===== Engineering Mathematics -- Homework 4, Question 11 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 20th March, 2026

# For the analysis  of equilibrium points near anchors

using DifferentialEquations
using LinearAlgebra
using Plots
using NLsolve

using Clustering # For NLsolve solution clustering
using ForwardDiff # To call the Jacobian of fixed points

# ---- Fundamental System Parameters ---- #

g = 9.810

# ---- Fundamental System Parameters ---- #

# A very small number to avoid 0 division
ϵ = 1e-6

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
    hat_A = vec_to_A / (dist_A + ϵ)
    hat_B = vec_to_B / (dist_B + ϵ)

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
    [2.0, 10.0],  # Anchor A
    [10.0, 10.0], # Anchor B
    0.0,          # Rest Length A
    0.0,          # Rest Length B
    10.0,         # k_A
    20.0,         # k_B
    0.0,          # c_A
    0.0,          # c_B
    0.5           # c_d (Air Drag)
)

# ---- Solving for a Static System ---- #

function static_forces!(F, r_c, sys::fparameters)
    a_vec = sys.rₐ - r_c
    b_vec = sys.rᵦ - r_c
    
    lₐ = norm(a_vec)
    lᵦ = norm(b_vec)
    
    hat_a = a_vec / lₐ
    hat_b = b_vec / lᵦ
    
    delta_a = lₐ - sys.Lₐ
    delta_b = lᵦ - sys.Lᵦ
    
    gravity = [0.0, -sys.m * g]
    force_a = sys.kₐ * delta_a * hat_a
    force_b = sys.kᵦ * delta_b * hat_b
    
    # Update the F array in-place
    F[1] = gravity[1] + force_a[1] + force_b[1] 
    F[2] = gravity[2] + force_a[2] + force_b[2] 
end

# ---- Generate Steady State Solution Clusters ---- #

num_samples = 2000
solutions_x = Float64[]
solutions_y = Float64[]

# Wrap the function so NLsolve only sees F and x
target_func!(F, x) = static_forces!(F, x, params)

for i in 1:num_samples
    x0 = -5.0 + 20.0 * rand()
    y0 = -10.0 + 25.0 * rand()
    initial_guess = [x0, y0]
    
    
    res = nlsolve(target_func!, initial_guess, autodiff=:forward, ftol=1.5)
    
    # If successful, save the raw, unrounded coordinates
    if converged(res)
        push!(solutions_x, res.zero[1])
        push!(solutions_y, res.zero[2])
    end
end

println("Generated $(length(solutions_x)) clustered equilibrium points.")
# ---- DO NOT Plot Results ---- #

# ---- Greedy Distance Grouping Algorithm ---- #
println("\n--- Starting Lightweight Distance Clustering ---")

points = tuple.(solutions_x, solutions_y)
clusters = Vector{Vector{Tuple{Float64, Float64}}}()

# Define how close points need to be to group together (meters)
tolerance = 1.5

for p in points
    found_cluster = false
    
    for cluster in clusters
        # Compare distance to the first point in the cluster (the seed)
        seed = cluster[1]
        dist = sqrt((p[1] - seed[1])^2 + (p[2] - seed[2])^2)
        
        if dist < tolerance
            push!(cluster, p)
            found_cluster = true
            break
        end
    end
    
    # If it didn't fit into any existing cluster, make a new one
    if !found_cluster
        push!(clusters, [p])
    end
end

println("Identified $(length(clusters)) unique equilibrium regions.")

# --- Calculate Centroids --- #
centroids_x = Float64[]
centroids_y = Float64[]

for (i, cluster) in enumerate(clusters)
    cx = sum(p[1] for p in cluster) / length(cluster)
    cy = sum(p[2] for p in cluster) / length(cluster)
    
    push!(centroids_x, cx)
    push!(centroids_y, cy)
    
    println("Region $i Centroid: (X: $(round(cx, digits=4)), Y: $(round(cy, digits=4)))")
end

println("Clustering automatically identified $(length(centroids_x)) physical equilibrium regions.")

final_solutions_plot = scatter(solutions_x, solutions_y, 
        label="Raw Iterations",
        size=(600, 500),
        markersize=1,
        markeralpha=0.3,
        markerstrokewidth=0,
        color=:blue,
        title="Final Steady State Regions & Centroids", 
        xlabel="X Coordinate", 
        ylabel="Y Coordinate",
        aspect_ratio=:equal
)

# Anchors for Context
scatter!(final_solutions_plot,
         [params.rₐ[1], params.rᵦ[1]], 
         [params.rₐ[2], params.rᵦ[2]], 
         label="Anchors A & B", 
         color=:red, 
         markersize=6, 
         markerstrokewidth=2
)

# Final Centroids
scatter!(final_solutions_plot,
         centroids_x, centroids_y,
         label="Steady State Centroids",
         color=:gold,
         markersize=8,
         markerstrokecolor=:black,
         markerstrokewidth=1.5,
         shape=:star
)

display(final_solutions_plot)

# ---- Analytic (almost) Analysis of Fixed Points ---- #

# Setup fixed points data structure
mutable struct fixedpoint_properties
    cx :: Float64
    cy :: Float64

    J_matrix :: Matrix{Float64}

    eigenvalues :: Vector{ComplexF64}

    stability :: String
end

# Setup array to hold fixed point data structures
fixed_points = Vector{fixedpoint_properties}()

# Wrapper Function for the system function
function get_derivatives(u)
    du = similar(u)

    spring_damper_mass_system!(du, u, params, 0.0)

    return du
end

# Populate fixed points array
for (cx, cy) in zip(centroids_x, centroids_y)
    fixed_point_state = [cx, cy, 0.0, 0.0]

    # Compute Jacobian
    J_matrix = ForwardDiff.jacobian(
        get_derivatives,            # System function
        fixed_point_state           # System state starting at the fixed point
    )

    #display(J_matrix)
    
    # Compute eigenvalues of the Jacobian matrix
    eigenvalues = eigvals(J_matrix)

    fixed_point = fixedpoint_properties(
        cx,
        cy,
        J_matrix,
        eigenvalues,
        ""
    )

    push!(
        fixed_points,
        fixed_point
    )
end

# -- Evaluate Stability of Collected Fixed Points -- #

# Zero tolerance
ztol = 1e-10

for i in eachindex(fixed_points)
    fixed_point = fixed_points[i]

    # Define flags
    negative_real = false
    positive_real = false
    zero_real = false

    has_imag = false
    zero_imag = false

    # Evaluate flags
    for λ in fixed_point.eigenvalues
        a = real(λ)
        b = imag(λ)

        if a > ztol
            positive_real = true
        elseif a < -ztol
            negative_real = true
        else
            zero_real = true
        end

        if abs(b) > ztol
            has_imag = true
        end     
    end

    # -- Evaluate based on Flags -- #
    # Evaluate based on the flags
    if negative_real && positive_real
        fixed_points[i].stability = "Saddle Point"

    elseif negative_real && !has_imag
        fixed_points[i].stability = "Stable Node"

    elseif negative_real && has_imag
        fixed_points[i].stability = "Stable Focus"
    
    elseif positive_real && !has_imag
        fixed_points[i].stability = "Unstable Node"

    elseif positive_real && has_imag
        fixed_points[i].stability = "Unstable Focus"
        
    elseif zero_real && !positive_real && !negative_real && has_imag
        fixed_points[i].stability = "Center (Marginally Stable)"
    
    else
        fixed_points[i].stability = "Unknown/Degenerate"
    
    end

end


# ---- Display Results (Pretty Print) ---- # Gemini 3.1

println("\n" * "="^50)
println(" SYSTEM EQUILIBRIUM ANALYSIS REPORT")
println("="^50)

for (i, fp) in enumerate(fixed_points)
    println("\nFixed Point [$i]:")
    println("  Coordinates : X = $(round(fp.cx, digits=4)), Y = $(round(fp.cy, digits=4))")
    println("  Classification : => $(uppercase(fp.stability)) <=")
    
    println("  Eigenvalues :")
    for (j, λ) in enumerate(fp.eigenvalues)
        r_part = round(real(λ), digits=4)
        i_part = round(imag(λ), digits=4)
        
        # Format nicely to handle the +/- im parts
        sign_str = i_part >= 0 ? "+" : "-"
        println("      λ$j = $r_part $sign_str $(abs(i_part))im")
    end
end
println("\n" * "="^50 * "\n")

# ---- Linearization about Fixed Points ---- #
# -- Define Common Parameters -- #
tspan_s = (0.0, 15.0)
tspan_u = (0.0,  2.0)
δu0 = [0.5, 0.5, 1.0, 2.0] # State Offset

# -- Linearization and Plot Block -- #
# ---- Linearization about Fixed Points ---- #

# Define linearized ODE function
function linearized_system!(du, u, p, t)
    J_matrix = p
    mul!(du, J_matrix, u)
end

# -- Modular Linearization and Plot Function -- #
function linearize_and_plot(
    fixed_point::fixedpoint_properties,
    type_label::String, δu0::Vector{Float64},
    tspan::Tuple{Float64, Float64}, 
    sys_params,
    tcolor::Symbol
    )
    
    # Extract point specific data
    J_point = fixed_point.J_matrix
    u_point = [fixed_point.cx, fixed_point.cy, 0.0, 0.0]

    # State setups
    u0_linear = δu0
    u0_nonlinear = δu0 + u_point

    # --- Simulate Linear --- #
    prob_lin = ODEProblem(linearized_system!, u0_linear, tspan, J_point)
    sol_lin = solve(prob_lin, Tsit5(), saveat=0.1)

    # --- Simulate Nonlinear --- #
    prob_nonlin = ODEProblem(spring_damper_mass_system!, u0_nonlinear, tspan, sys_params)
    sol_nonlin = solve(prob_nonlin, Tsit5(), saveat=0.1)

    # ---- Visual Proof (Overlay Plot) ---- # by Gemini 3.1 Pro
    
    # 1. Extract Nonlinear data
    nonlin_x = [u[1] for u in sol_nonlin.u]
    nonlin_y = [u[2] for u in sol_nonlin.u]

    # 2. Extract Linear data and shift it
    lin_x = [u[1] + u_point[1] for u in sol_lin.u]
    lin_y = [u[2] + u_point[2] for u in sol_lin.u]

    # 3. Create the base plot
    comparison_plot = plot(
        nonlin_x, nonlin_y, 
        label="Nonlinear (True Physics)", 
        lw=2, color=tcolor, 
        title="Local Dynamics: $(type_label) Point",
        xlabel="X Position (m)", ylabel="Y Position (m)",
        aspect_ratio=:equal
    )

    # 4. Overlay the shifted Linear path
    plot!(comparison_plot, 
        lin_x, lin_y, 
        label="Linearized (Jacobian Approx)", 
        lw=2, linestyle=:dash, color=:orange
    )

    # 5. Mark the Fixed Point
    scatter!(comparison_plot, 
        [u_point[1]], [u_point[2]], 
        label="$(type_label) Fixed Point", 
        shape=:star, color=:gold, 
        markersize=8, markerstrokewidth=1
    )

    display(comparison_plot)
end

# -- Linearization about a stable fixed point -- #

stable_point_index = -1

for i in eachindex(fixed_points)
    stability = fixed_points[i].stability
    if stability == "Stable Node" || stability == "Stable Focus"
        global stable_point_index = i
        break
    end
end

println("Generating comparison plot for Stable Point at index $stable_point_index...")

stable_point = fixed_points[stable_point_index]
#
linearize_and_plot(
    stable_point,       # The struct
    "Stable",           # The label string
    δu0,                # Offset
    tspan_s,            # Time
    params,             # The system parameters
    :blue               # Line colour
)
#

# -- Linearization about an unstable fixed point -- #
#=
unstable_point_index = -1

for i in eachindex(fixed_points)
    stability = fixed_points[i].stability
    if stability != "Stable Node" && stability != "Stable Focus"
        global unstable_point_index = i
        break
    end
end

println("Generating comparison plot for Unstable Point at index $unstable_point_index...")


unstable_point = fixed_points[unstable_point_index]

linearize_and_plot(
    unstable_point,       # The struct
    "Unstable",           # The label string
    δu0 + ,                  # Offset
    tspan_u,              # Time
    params,               # The system parameters
    :red                  # Line colour
)
=#
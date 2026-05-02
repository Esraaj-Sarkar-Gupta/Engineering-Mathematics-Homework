# ===== Engineering Mathematics -- Homework 3, Question 9 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 15th February, 2026

using DifferentialEquations
using LinearAlgebra
using Plots
using NLsolve

# ---- Fundamental System Parameters ---- #

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
    [2.0, 10.0],  # Anchor A
    [10.0, 10.0], # Anchor B
    0.0,          # Rest Length A
    0.0,          # Rest Length B
    10.0,         # k_A
    20.0,         # k_B
    10.0,         # c_A
    1.0,          # c_B
    0.5           # c_d (Air Drag)
)
# Initial Conditions
u0 = [7.3333, 8.3649, 1.0, 2.0] 

# Time Span
tspan = (0.0, 50.0)

prob = ODEProblem(spring_damper_mass_system!, u0, tspan, params)
sol = solve(prob, Tsit5(), saveat=0.1)

phase_space = plot(sol, idxs=(1,2), 
     title="Mass Trajectory", 
     xlabel="x position (m)", 
     ylabel="y position (m)", 
     label="Path", lw=2, aspect_ratio=:equal)

# Add anchors to the plot for context
scatter!(phase_space, [params.rₐ[1], params.rᵦ[1]], [params.rₐ[2], params.rᵦ[2]], 
         label="Anchors", markersize=8, color=:red)

display(phase_space)

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
# ---- Plot Results ---- #

# Add the anchors to the plot
static_solutions = scatter(
         [params.rₐ[1], params.rᵦ[1]], 
         [params.rₐ[2], params.rᵦ[2]], 
         label="Anchors A & B", 
         color=:red, 
         markersize=6, 
         markerstrokewidth=2)

scatter!(static_solutions,
        solutions_x, solutions_y, 
        label="Equilibrium Points",
        size=(500,500),
        markersize=3,
        color=:blue,
        title="Static Equilibrium Positions", 
        xlabel="X Coordinate", 
        ylabel="Y Coordinate",
        aspect_ratio=:equal)

display(static_solutions)

# ---- Custom Steady State Clustering Algorithm ---- #

# Combine points for easier handling
points = tuple.(solutions_x, solutions_y)

# Round off to integers and find unique classes
rounded_points = tuple.(round.(Int, solutions_x), round.(Int, solutions_y))
unique_classes = unique(rounded_points)

println("Identified $(length(unique_classes)) integer class seeds.")

# Create a dictionary to hold the assigned points for each class
clusters = Dict(c => Tuple{Float64, Float64}[] for c in unique_classes)

# Assign each original point to the nearest integer class using Euclidean distance
for p in points
    best_class = unique_classes[1]
    min_dist = Inf
    
    for c in unique_classes
        # Euclidean distance formula
        dist = sqrt((p[1] - c[1])^2 + (p[2] - c[2])^2)
        if dist < min_dist
            min_dist = dist
            best_class = c
        end
    end
    push!(clusters[best_class], p)
end

# Centroids of each newly formed class
centroids_x = Float64[]
centroids_y = Float64[]

println("\n--- Final Steady State Regions ---")
for (class_seed, cluster_points) in clusters
    if !isempty(cluster_points)
        # Calculate mean of X and Y
        cx = sum(p[1] for p in cluster_points) / length(cluster_points)
        cy = sum(p[2] for p in cluster_points) / length(cluster_points)
        
        push!(centroids_x, cx)
        push!(centroids_y, cy)
        
        println("Region Centroid: (X: $(round(cx, digits=4)), Y: $(round(cy, digits=4)))")
    end
end

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

# ---- Stability Check via Dynamic Simulation ---- #

println("\n--- Running Dynamic Stability Check ---")

# Create new plot to hold the stability of discovered fixed points

# Add anchors and centroids for visual reference
stability_plot = scatter([params.rₐ[1], params.rᵦ[1]], [params.rₐ[2], params.rᵦ[2]], 
         label="Anchors A & B", color=:red, markersize=6)

plot!(stability_plot, title="Stability Check at Centroids", 
                      xlabel="X Coordinate (m)", 
                      ylabel="Y Coordinate (m)", 
                      aspect_ratio=:equal, 
                      size=(600, 500))
         
scatter!(stability_plot, centroids_x, centroids_y, 
         label="Centroids", color=:gold, shape=:star, markersize=8, markerstrokewidth=1)

# Define velocity perturbations
v_perturb_x = 0.1
v_perturb_y = 0.1

# Loop through each discovered centroid
for i in 1:length(centroids_x)
    u0_stab = [centroids_x[i], centroids_y[i], v_perturb_x, v_perturb_y]
    
    tspan_stab = (0.0, 20.0) 
    
    # Solve the ODE
    prob_stab = ODEProblem(spring_damper_mass_system!, u0_stab, tspan_stab, params)
    sol_stab = solve(prob_stab, Tsit5(), saveat=0.1)
    
    # Plot the resulting path
    plot!(stability_plot, sol_stab, idxs=(1,2), 
          label="Trajectory $i", lw=2, alpha=0.8)
end

# Display the final stability plot
display(stability_plot)
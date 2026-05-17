# ===== Engineering Mathematics -- Question 21 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 2nd May, 2026

# Optimal Control: Deterministic Backward Reachability Boundary

using Plots

g = 1.62

println("---- Run ---- ")

# ---- Deterministic Boundary Mapping ---- #
function map_reachability_boundary(
    num_angles, 
    T, 
    N, 
    u_max
)
    h = T / N
    
    # Arrays to store the discovered t=0 boundary states
    boundary_x0 = Float64[]
    boundary_z0 = Float64[]

    println("Starting Deterministic Boundary Sweep ($num_angles angles)...")

    # Sweep theta from 0 to 2π
    angles = range(0, 2*pi, length=num_angles)

    for (i, theta) in enumerate(angles)
        # Start at the target origin at T_final
        x, z, vx, vz = 0.0, 0.0, 0.0, 0.0
        
        valid_trajectory = true

        # Lock the max control effort in the current direction theta
        ux = u_max * cos(theta)
        uz = u_max * sin(theta)

        # Step backward in time
        for k in 1:N
            # Reverse Euler Integration
            vx_prev = vx - h * ux
            vz_prev = vz - h * (uz - g)
            
            x_prev = x - h * vx_prev
            z_prev = z - h * vz_prev

            # Update state for the next backward step
            x = x_prev
            z = z_prev
            vx = vx_prev
            vz = vz_prev

            # State constraint: Check if trajectory goes underground
            if z < 0.0
                valid_trajectory = false
                break
            end
        end

        # Only save the start points that never violate the ground constraint
        if valid_trajectory
            push!(boundary_x0, x)
            push!(boundary_z0, z)
        end
        
        if i % 100 == 0
            println("Swept $i / $num_angles directions...")
        end
    end

    println("Finished! Mapped $(length(boundary_x0)) boundary points.")
    return boundary_x0, boundary_z0
end

# ---- Parameter Setup ---- #
T_final = 5.0
N_steps = 100
u_max_limit = 100.0
sweep_resolution = 1000 # Number of angles to test

# Extract the boundary points
x0_bounds, z0_bounds = map_reachability_boundary(sweep_resolution, T_final, N_steps, u_max_limit)

# ---- Plotting the Boundary ---- #
# We use a line plot (`plot`) instead of a scatter plot because we are drawing a continuous perimeter
boundary_plot = plot(
    x0_bounds, z0_bounds, 
    label="Max Effort Boundary", 
    color=:blue, 
    linewidth=2,
    fill=(0, 0.2, :blue), # Fills the area inside the boundary with a light blue
    title="Backward Reachable Set (Constant Control Angle)", 
    xlabel="Initial X Position [m]", 
    ylabel="Initial Z Position [m]",
    legend=:outerright,
    size=(800, 600),
    aspect_ratio=:equal # Keeps the geometry proportional
)

# Mark the target destination
scatter!(boundary_plot, [0.0], [0.0], color=:red, markersize=5, label="Target Destination (0,0)")

display(boundary_plot)
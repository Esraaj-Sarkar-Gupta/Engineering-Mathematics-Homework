# ===== Engineering Mathematics -- Question 21 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 2nd May, 2026

# Optimal Control: Backward Reachability Analysis

using Plots
using Distributions

using Statistics
using KernelDensity

g = 1.62

println("---- Run ---- ")

# ---- Backward Integration Function ---- #
function backward_random_reachability(
    num_samples, 
    T, 
    N, 
    u_max
)
    h = T / N
    
    # Arrays to store the discovered t=0 starting states
    valid_x0 = Float64[]
    valid_z0 = Float64[]
    valid_vx0 = Float64[]
    valid_vz0 = Float64[]

    println("Starting Backward Reachability Sampling ($num_samples trajectories)...")

    for i in 1:num_samples
        # Start at the end.
        #  Life Can Only Be Understood Backwards, But It Must Be Lived Forwards.
        x = 0.0
        z = 0.0
        vx = 0.0
        vz = 0.0
        
        valid_trajectory = true

        # Step backward in time
        for k in 1:N
            theta = rand(Uniform(0, 2*pi))
            #theta = k / (2*pi*N)
            r = u_max * sqrt(rand()) 
            
            ux = r * cos(theta)
            uz = r * sin(theta)

            # Reverse Euler Integration
            # We must calculate v_prev first to use it for x_prev, exactly mirroring your forward loop
            vx_prev = vx - h * ux
            vz_prev = vz - h * (uz - g)
            
            x_prev = x - h * vx_prev
            z_prev = z - h * vz_prev

            # Update state for the next backward step
            x = x_prev
            z = z_prev
            vx = vx_prev
            vz = vz_prev

            # State constraint: The bot cannot travel underground at any point in the timeline
            if z < 0.0
                valid_trajectory = false
                break
            end
        end

        if valid_trajectory
            push!(valid_x0, x)
            push!(valid_z0, z)
            push!(valid_vx0, vx)
            push!(valid_vz0, vz)
        end
        
        # Progress tracker
        if i % 250000 == 0
            println("Simulated $i / $num_samples backward trajectories...")
        end
    end

    println("Finished! Found $(length(valid_x0)) viable starting points.")
    return valid_x0, valid_z0, valid_vx0, valid_vz0
end

# ---- Running it all ---- #
T_final = 5.0
N_steps = 100
u_max_limit = 10.0
samples = 1e8 # Cheap computation

# Extract the parameter space of valid t=0 starting conditions
x0, z0, vx0, vz0 = backward_random_reachability(samples, T_final, N_steps, u_max_limit)

# ---- Plotting the Boundary ---- #
scatter_plot = scatter(
    x0, z0, 
    label="Viable Starts (t=0)", 
    color=:purple, 
    markersize=4, 
    markerstrokewidth=0,
    alpha=0.8,
    title="Backward Reachability: Viable Initial Positions", 
    xlabel="Initial X Position [m]", 
    ylabel="Initial Z Position [m]",
    legend=:outerright,
    size=(800, 600)
)

display(scatter_plot)

println("Calculating statistics and KDE for $(length(x0)) points...")

mean_x = mean(x0)
mean_z = mean(z0)
std_x = std(x0)
std_z = std(z0)

println("Mean Initial X: ", round(mean_x, digits=2), " ± ", round(std_x, digits=2), " m")
println("Mean Initial Z: ", round(mean_z, digits=2), " ± ", round(std_z, digits=2), " m")

kde_result = kde((x0, z0))

contour_plot = contour(
    kde_result.x, kde_result.y, kde_result.density,
    fill=true, 
    c=:viridis,
    levels=15,
    title="Backward Reachability: Interpolated Density",
    xlabel="Initial X Position [m]",
    ylabel="Initial Z Position [m]",
    colorbar_title="Point Density",
    size=(800, 600)
)

scatter!(contour_plot, [mean_x], [mean_z], 
    color=:red, markersize=6, marker=:star, label="Mean Start"
)

plot!(contour_plot, 
    [mean_x - std_x, mean_x + std_x, mean_x + std_x, mean_x - std_x, mean_x - std_x], # X coordinates for the box
    [mean_z - std_z, mean_z - std_z, mean_z + std_z, mean_z + std_z, mean_z - std_z], # Z coordinates for the box
    color=:white, linestyle=:dash, lw=2, label="1-Sigma Deviation"
)

display(contour_plot)
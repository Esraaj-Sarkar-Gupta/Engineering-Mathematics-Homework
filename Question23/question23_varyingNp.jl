# ===== Engineering Mathematics -- Question 23 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 17th May, 2026
# Open-Loop Collocation tied into Model Predictive Control

using JuMP
using Ipopt
using Plots
using LinearAlgebra

# ---- Global System Parameters ---- #
tspan = (0.0, 5.0)    
step_t = 10 * 1e-3  # 10 ms
n_steps = Int((tspan[2] - tspan[1]) / step_t)

state_initial = [5.0, 10.0, -1.0, -2.0]
state_end     = [0.0, 0.0, 0.0, 0.0]
g = 1.62
u_max = 10.0

# MPC Tuning Parameters
Np = 10               # Prediction Horizon
Q = Diagonal([10.0, 10.0, 1.0, 1.0]) 
R = Diagonal([0.1, 0.1])

# Gaussian Distribution for Noise
μ = 0.80
σ = 0.33



# ---- Reference Generator (Question 21) ---- #
println("Generating Open-Loop Reference Trajectory....")

function optimize_trajectory(N, s0, sT, tspan, u_max, h)
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    @variables(model, begin
        Xx[1:N+1]; Xz[1:N+1]
        vx[1:N+1]; vz[1:N+1]
        ux[1:N+1]; uz[1:N+1]
    end)

    # Trapezoid Collocation
    for k in 1:N
        @constraint(model, Xx[k+1] == Xx[k] + (h/2) * (vx[k] + vx[k+1]))
        @constraint(model, Xz[k+1] == Xz[k] + (h/2) * (vz[k] + vz[k+1]))
        @constraint(model, vx[k+1] == vx[k] + (h/2) * (ux[k] + ux[k+1]))
        @constraint(model, vz[k+1] == vz[k] + (h/2) * (uz[k] + uz[k+1] - 2*g))
    end

    @objective(model, Min, sum(sqrt(ux[k]^2 + uz[k]^2 + 1e-6) * h for k in 1:N+1))

    @constraint(model, Xx[1] == s0[1]);  @constraint(model, Xx[N+1] == sT[1])
    @constraint(model, Xz[1] == s0[2]);  @constraint(model, Xz[N+1] == sT[2])
    @constraint(model, vx[1] == s0[3]);  @constraint(model, vx[N+1] == sT[3])
    @constraint(model, vz[1] == s0[4]);  @constraint(model, vz[N+1] == sT[4])

    @constraint(model, Xz .>= 0)
    @constraint(model, [k=1:N+1], ux[k]^2 + uz[k]^2 <= u_max^2)

    optimize!(model)
    return value.(Xx), value.(Xz), value.(vx), value.(vz)
end

# Run the open-loop solver to get the ideal arrays
xx_ref, xz_ref, vx_ref, vz_ref = optimize_trajectory(n_steps, state_initial, state_end, tspan, u_max, step_t)

# Package the arrays into a vector of state vectors for easy MPC horizon slicing
ref_trajectory = [[xx_ref[k], xz_ref[k], vx_ref[k], vz_ref[k]] for k in 1:n_steps+1]

# ---- MCP Controller ---- #
println("Running Closed-Loop MPC Tracker....")

function solve_mcp(state_current, state_ref_horizon, Q, R, Np, h, u_max)
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    @variables(model, begin
        x[1:Np+1]; z[1:Np+1]
        vx[1:Np+1]; vz[1:Np+1]
        ux[1:Np+1]; uz[1:Np+1]
    end)
    
    @constraints(model, begin
        x[1] == state_current[1]
        z[1] == state_current[2]
        vx[1] == state_current[3]
        vz[1] == state_current[4]
    end)

    for k in 1:Np
        @constraint(model, x[k+1] == x[k] + (h/2)*(vx[k] + vx[k+1]))
        @constraint(model, z[k+1] == z[k] + (h/2)*(vz[k] + vz[k+1]))
        @constraint(model, vx[k+1] == vx[k] + (h/2)*(ux[k] + ux[k+1]))
        @constraint(model, vz[k+1] == vz[k] + (h/2)*(uz[k] + uz[k+1] - 2*g))
    end

    @constraint(model, [k = 1:Np+1], z[k] >= 0)
    @constraint(model, [k = 1:Np+1], ux[k]^2 + uz[k]^2 <= u_max^2)

    # -- Define Objective -- #
    q = Q.diag
    r = R.diag
    
    @objective(model, Min,
        sum(
            q[1] * (x[k]  - state_ref_horizon[k][1])^2 +
            q[2] * (z[k]  - state_ref_horizon[k][2])^2 +
            q[3] * (vx[k] - state_ref_horizon[k][3])^2 +
            q[4] * (vz[k] - state_ref_horizon[k][4])^2 +
            r[1] * ux[k]^2 +
            r[2] * uz[k]^2
            for k in 1:Np+1
        )
    )

    optimize!(model)
    
    if termination_status(model) ∉ (MOI.LOCALLY_SOLVED, MOI.OPTIMAL) # Julia is so cool
        return 0.0, 0.0
    end

    return value(ux[1]), value(uz[1])
end

# Plant step function
function step_dynamics(s, u, h, g)
    # Unpack state vectors
    px, pz, vx, vz = s
    ux, uz = u

    # Using simple Euler for the forward step
    return [px + h*vx, pz + h*vz, vx + h*ux, vz + h*(uz - g)]
end

# ---- Main Execution Loop (Function) ---- #

function simulate_mpc_run(Np, n_steps, state_initial, ref_trajectory, Q, R, step_t, u_max, g, μ, σ)
    println("Running MPC Simulation for Np = $Np...")
    state_actual = copy(state_initial)
    trajectory_mpc = [copy(state_actual)]
    
    for t in 1:n_steps
        reference_start = t
        reference_end = min(t + Np, n_steps + 1)
        
        # Slice the dynamically feasible path
        state_ref_horizon = ref_trajectory[reference_start:reference_end]

        # Pad with the final target state
        while length(state_ref_horizon) < Np + 1
            push!(state_ref_horizon, ref_trajectory[end])
        end

        # Calculate optimal control
        ux_opt, uz_opt = solve_mcp(state_actual, state_ref_horizon, Q, R, Np, step_t, u_max)

        # Inject the Gaussian disturbance
        α = μ + σ * randn()
        α = clamp(α, 0.0, 1.0) # Keep efficiency strictly between

        u_applied = [α * ux_opt, α * uz_opt]

        # Step the true physics forward with the degraded control
        state_actual = step_dynamics(state_actual, u_applied, step_t, g)
        push!(trajectory_mpc, copy(state_actual))

        # Check for crashes
        if state_actual[2] <= 0.0
            break
        end
    end
    
    return trajectory_mpc
end

# ---- Execute and Plot ---- #

# Horizon lengths to test (from your homework specs)
Np_list = [5, 10, 20, 30]

# Setup the base plot (Large size as requested)
p1 = plot(xx_ref, xz_ref, 
          label="Open-Loop Ref", 
          lw=2, 
          ls=:dash, 
          color=:black, 
          size=(1000, 800), # Large graph
          legend=:bottomleft)

# Define distinct colors for the different horizons
colors = [:crimson, :steelblue, :darkorange, :purple]

# Run the simulation for each Np and add it to the plot
for (i, current_Np) in enumerate(Np_list)
    traj = simulate_mpc_run(current_Np, n_steps, state_initial, ref_trajectory, Q, R, step_t, u_max, g, μ, σ)
    
    # Extract coordinates
    x_m = [s[1] for s in traj]
    z_m = [s[2] for s in traj]
    
    # Plot with alpha = 0.8
    plot!(p1, x_m, z_m, 
          label="MPC (Np = $current_Np)", 
          lw=2, 
          color=colors[i], 
          alpha=0.8)
    
    # Mark the end point (whether target reached or crashed)
    scatter!(p1, [x_m[end]], [z_m[end]], label="", color=colors[i], ms=5)
end

# Add start marker
scatter!(p1, [state_initial[1]], [state_initial[2]], label="Start", color=:green, ms=8)

xlabel!("x [m]")
ylabel!("z [m]")
title!("Trajectory Tracking vs. Prediction Horizon (Np)")

display(p1)
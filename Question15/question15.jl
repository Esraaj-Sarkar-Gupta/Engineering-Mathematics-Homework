# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 23rd April, 2026

using Plots

# ---- System Definition ---- #
function system!(u)
    x, v = u
    dv = -1
    dx = v
    return dx, dv
end

# ---- Integration Methods ---- #
function euler_integration(u, time_step = 0.01)
    x, v = u
    dx, dv = system!(u)
    v_new = v + time_step * dv
    x_new = x + time_step * dx
    return x_new, v_new
end

function trapezoid_integration(u, time_step = 0.01)
    x, v = u
    dx, dv = system!(u)
    v_new = v + time_step * dv 
    x_new = x + (time_step / 2.0) * (v + v_new)
    return x_new, v_new
end

# ---- Analytical Solution ---- #
exact_v(t) = 0.5 - t
exact_x(t) = 0.5 * t - 0.5 * t^2

# ---- Simulation Function ---- #
# Wrapping the loop in a function avoids global scope issues and makes it reusable
function run_simulation(h, time_end=20)
    iterations = Int(time_end / h)
    
    state_euler = (0.0, 0.5)
    error_euler = Float64[]
    
    state_trap = (0.0, 0.5)
    error_trap = Float64[]
    
    time_array = Float64[]
    
    for iteration in 1:iterations
        current_time = iteration * h
        push!(time_array, current_time)
        
        x_true = exact_x(current_time)
        v_true = exact_v(current_time)
        
        # -- Euler Step and Error -- #
        state_euler = euler_integration(state_euler, h)
        e_x_euler = state_euler[1] - x_true
        e_v_euler = state_euler[2] - v_true
        rms_euler = sqrt((e_x_euler^2 + e_v_euler^2) / 2.0)
        push!(error_euler, rms_euler)
        
        # -- Trapezoid Step and Error -- #
        state_trap = trapezoid_integration(state_trap, h)
        e_x_trap = state_trap[1] - x_true
        e_v_trap = state_trap[2] - v_true
        rms_trap = sqrt((e_x_trap^2 + e_v_trap^2) / 2.0)
        push!(error_trap, rms_trap)
    end

    # Create the Comparison Plot
    p_comp = plot(time_array, error_euler, label = "Euler", color = "red", linewidth = 2,
                  title = "Euler vs Trap (h = $h)", xlabel = "Time (s)", ylabel = "RMS Error",
                  grid = true, legend = :topleft)
    plot!(p_comp, time_array, error_trap, label = "Trapezoid", color = "blue", linewidth = 2)

    # Create the Trapezoid-Only Plot
    p_trap = plot(time_array, error_trap, label = "Trapezoid", color = "blue", linewidth = 2,
                  title = "Trap Error Only (h = $h)", xlabel = "Time (s)", ylabel = "RMS Error",
                  grid = true, legend = :topleft)
                  
    return p_comp, p_trap
end

# ---- Run Simulations for different 'h' values ---- #
h_values = [1.0, 0.1]

# Run for h = 1.0
p_comp_1, p_trap_1 = run_simulation(0.01)

# Run for h = 0.1
p_comp_01, p_trap_01 = run_simulation(0.001)

# ---- Combine into a 2x2 Grid ---- #
final_plot = plot(
    p_comp_1, p_comp_01, 
    p_trap_1, p_trap_01, 
    layout = (2, 2), 
    size = (1000, 800),
    margin = 5Plots.mm
)

display(final_plot)

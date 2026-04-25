# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 25th April, 2026

# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 23rd April, 2026

# Solving Boundary Value Problems (BVPs)

using Plots
using Printf

# ---- Define the System ---- #
x_start = 0.0
x_end = 0.0 

t_start = 0.0
t_end = 1.0

function system!(u)
    x, v = u

    dv = -1
    dx = v

    du = dx, dv
    return du
end

# ---- Define Solvers ---- #
# Borrowed from Question 15

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

# ---- System Parameters ---- #
# Temporal
time_step = 0.01
end_time = 1
iterations = end_time / time_step

# Boundary Conditions
x0 = 0
xend = 0

function shooting_function(s)
    v0 = s

    state = x0, v0
    for iteration in 1:iterations
        state = euler_integration(state, time_step)
    end

    return state[1] - xend # Return F(s) = x(end_time; s) - x(T)
end

# ---- Define Bijection Root-Finder ---- #
function bijection_solve(
    func,
    a,
    b,
    tol = 1e-5,
    max_iter = 1e+5
)
    fa = func(a)
    fb = func(b)

    # NOTE: f(a) and f(b) must have different signs in order for the method to work
    # Check for the above Conditions

    if (fa * fb > 0)
        return true  # Return true :: Error flag = 1
    end

    # Init midpoint
    c = a
    
    for i in 1:max_iter
        c = (a+b)/2
        fc = func(c)

        if (abs(fc) < tol) || ((b - a) / 2 < tol)
            return c
        end

        # Check for sign changes
        if (fa * fc < 0)
            # Discard right hand side
            fb = fc 
            b = c
        elseif (fb * fc < 0)
            # Discard left hand side
            fa = fc 
            a = c
        end
    end

    return c
end

# ---- Finding s ---- #
# Wrap the shooting function with the root Finder

result = bijection_solve(
    shooting_function,
    -10,
    10,
    1e-5,
    10000
)

if (result === true)
    println("Encountered Error!")
else
    @printf("v(0) = %f", result)
end

# ---- Simulate System using Result ---- #
state = [x_start, result]

h = 0.01
iterations = Int((t_end - t_start) / h)

# Arrays to store the shooting trajectory
t_shooting = range(t_start, t_end, length = iterations + 1)
x_shooting = zeros(iterations + 1)
v_shooting = zeros(iterations + 1)

# Initialize start values
x_shooting[1] = state[1]
v_shooting[1] = state[2]

for i in 1:iterations
    # Using Euler integration for uniformity as requested
    global state = euler_integration(state, h)
    x_shooting[i+1] = state[1]
    v_shooting[i+1] = state[2]
end

# ---- Collocation Method ---- #

function collocation_solve(N, t_start = 0.0, t_end = 1.0)
    h_colloc = (t_end - t_start) / N
    num_points = N + 1

    # Extract system acceleration (dx goes to _, dv goes to a)
    _, a = system!((0.0, 0.0))

    # Total unknowns
    total_vars = 2 * num_points

    # Init matrices
    A = zeros(total_vars, total_vars)
    b = zeros(total_vars) 

    # x(0) = 0
    A[1,1] = 1.0
    b[1] = 0.0

    # x(N+1) = 0
    A[2, num_points] = 1.0
    b[2] = 0.0

    # -- Enforce Euler Discretization -- #
    row = 3 # 1 and 2 have been used up in the boundary Conditions

    for k in 1:N
        # -- Index Helpers -- #
        idx_xk = k
        idx_xnext = k + 1
        idx_vk = num_points + k
        idx_vnext = num_points + k + 1

        # -- Populate Matrices -- #

        # Rule: x_{k+1} - x_k - h*v_k = 0
        A[row, idx_xnext] = 1.0
        A[row, idx_xk] = -1.0
        A[row, idx_vk] = -h_colloc
        b[row] = 0.0

        row += 1

        # Rule: v_{k+1} - v_k = a*h
        A[row, idx_vnext] = 1.0
        A[row, idx_vk] = -1.0
        b[row] = a * h_colloc

        row += 1
    end

    # --- Solve for the Linear System Az = b ---- #
    z = A \ b

    # ---- Unpack Solution ---- #
    x_sol = z[1:num_points]
    v_sol = z[num_points+1:end]
    time_array = range(t_start, t_end, length = num_points)

    return time_array, x_sol, v_sol
end

# ---- Run Collocation ---- #
# Match the step size to the shooting method
N_colloc = iterations 
t_colloc, x_colloc, v_colloc = collocation_solve(N_colloc, t_start, t_end)

# ---- Plot Comparison ---- #
p_compare = plot(
    t_shooting, x_shooting,
    label = "Shooting Method",
    linewidth = 5,
    color = :green,
    alpha = 0.5,
    title = "Shooting vs Collocation Comparison",
    xlabel = "t",
    ylabel = "x(t)",
    grid = true,
    legend = :topleft
)

plot!(
    p_compare,
    t_colloc, x_colloc,
    label = "Collocation Method",
    linewidth = 5,
    color = :orange,
    alpha = 0.5,
    linestyle = :dash # Added a dash so you can actually see both overlapping lines
)

display(p_compare)

# ---- Compute Residuals ---- #
residuals = x_shooting .- x_colloc

# ---- Plot Residuals ---- #
p_residual = plot(
    t_shooting, residuals,
    label = "x_shot - x_colloc",
    linewidth = 2,
    color = :purple,
    title = "Numerical Residual Over Time",
    xlabel = "t",
    ylabel = "Residual (m)",
    grid = true,
    legend = :topright
)

# Adding a horizontal line at 0 for reference
hline!([0], color = :black, linestyle = :dash, label = "")

display(p_residual)

# Just for your own sanity check in the console
max_err = maximum(abs.(residuals))
@printf("\nMaximum residual between methods: %.2e m\n", max_err)






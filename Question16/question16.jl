# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 23rd April, 2026

# Solving Boundary Value Problems (BVPs)

using Plots
using Printf

# ---- Define the System ---- #
function system!(u)
    x, v = u

    dv = -1
    dx = v

    du = dx, dv
    return du
end

# ---- Define Solver ---- #
# Borrowed from Question 15
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
        state = trapezoid_integration(state, time_step)
    end

    return state[1] - xend # Return F(s) = x(end_time; s) - x(T)
end

# ---- Run Iterations of F vs s ---- #
s_values = range(-5, 5, 1000)
F_values :: Array = []

for s in s_values
    F = shooting_function(s)
    
    push!(F_values, F)
end

F_vs_s_plot = plot(
    s_values,
    F_values,
    xlabel = "s",
    ylabel = "F(s)",
    title = "F vs s Plot",
    legend = false,
    grid = true,
    color = :green
)

scatter!(
    F_vs_s_plot,
    [1],
    [0],
    color =:red
)

display(F_vs_s_plot)

# ---- Define bisection Root-Finder ---- #
function bisection_solve(
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

result = bisection_solve(
    shooting_function,
    -100,
    100,
    1e-5,
    100000
)

if (result === true)
    println("Encountered Error!")
else
    @printf("v(0) = %f", result)
end








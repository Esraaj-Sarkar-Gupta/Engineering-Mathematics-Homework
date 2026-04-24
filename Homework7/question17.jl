# ===== Engineering Mathematics -- Homework 7 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 25th April, 2026

# Collocation

using Plots

function system!(u)
    x, v = u

    a = -1

    dv = a
    dx = v

    du = dx, dv
    return du, a
end

function euler_integration(u, time_step = 0.01)
    x, v = u
    dx, dv = system!(u)
    v_new = v + time_step * dv
    x_new = x + time_step * dx
    return x_new, v_new
end

# ---- Collocation Method ---- #

function collocation_solve(N, t_start = 0.0, t_end = 1.0)
    h = (t_end - t_start) / N
    num_points = N + 1

    # Extract system acceleration
    _, a = system!((0,0))

    # Total unknowns
    total_vars = 2 * num_points

    # Init matrices
    A = zeros(total_vars, total_vars)
    b = zeros(total_vars) # Collocation

    # x(1) = 0
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
        A[row, idx_vk] = -h
        b[row] = 0.0

        row += 1

        # Rule: v_{k+1} - v_k = a*h
        A[row, idx_vnext] = 1.0
        A[row, idx_vk] = -1.0
        b[row] = a*h

        row += 1
    end

    # --- Solve for the Linear System Az = b ---- #
    z = A \ b

    # ---- Upack Solution ---- #
    x_sol = z[1:num_points]
    v_sol = z[num_points+1:end]
    time_array = range(t_start, t_end, length = num_points)

    return time_array, x_sol, v_sol
end

# ---- Run and Plot Results ---- #
N = 100
t, x, v = collocation_solve(N)

p_colloc_xt = plot(
    t, x,
    labels = "x(t)",
    #marker = :cross,
    title = "Collocation Method: N = $N",
    xlabel = "t",
    ylabel = "x(t)",
    grid = true,
    legend = :topleft,
    color = :green
)

p_colloc_xv = plot(
    x, v,
    #labels = "",
    #marker = :cross,
    title = "Collocation Method: N = $N",
    xlabel = "x(t)",
    ylabel = "v(t)",
    grid = true,
    legend = false,
    color = :orange
)

display(p_colloc_xt); display(p_colloc_xv)
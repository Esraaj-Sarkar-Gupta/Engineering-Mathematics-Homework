# ===== Engineering Mathematics -- Homework 4, Question 10 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 20th March, 2026

using DifferentialEquations
using LinearAlgebra
using Plots

# --- Define System ---- #
A_mat = [
    3.0 1.0;
    1.0 3.0
]
X0 = [5.0, 0.0]
tspan = (0.0, 2.0) 

# ---- Analytical Solution ---- #
function analytical_solution(t)
    v1 = [1.0, 1.0]
    v2 = [1.0, -1.0]
    λ1, λ2 = 4.0, 2.0
    c1 = c2 = 2.5

    # FIXED: Added explicit multiplication * before vectors
    return c1 * exp(λ1 * t) * v1 + c2 * exp(λ2 * t) * v2
end

# ---- Numerical Solution ---- #
function system!(dX, X, p, t)
    # p is the matrix A passed from ODEProblem
    mul!(dX, p, X)
end

# FIXED: Corrected argument order (f, u0, tspan, p)
problem = ODEProblem(system!, X0, tspan, A_mat)
solution = solve(problem, Tsit5(), saveat=0.05)

# -- Extract Numerical Solutions -- #
t_vals = solution.t
z1_num = [u[1] for u in solution.u]
z2_num = [u[2] for u in solution.u]

# -- Generate Analytical Solutions -- #
mutable struct anal_sol
    t :: Vector{Float64}
    u :: Matrix{Float64} # Changed to Matrix for easier indexing
end

# Initialize the struct with the time vector and a correctly sized matrix
n_steps = length(t_vals)
u_anal = zeros(n_steps, 2)

for i in 1:n_steps
    res = analytical_solution(t_vals[i])
    u_anal[i, 1] = res[1]
    u_anal[i, 2] = res[2]
end

# Instantiate your struct
comparison_data = anal_sol(t_vals, u_anal)

# ---- Plotting Results ---- #
p1 = plot(title="Analytical vs Numerical: ż = Az", xlabel="Time (s)", ylabel="State")

# Plot Z1 comparison
plot!(p1, t_vals, z1_num, label="Numerical z1", lw=4, color=:lightblue)
plot!(p1, comparison_data.t, comparison_data.u[:, 1], 
      label="Analytical z1", lw=2, linestyle=:dash, color=:blue)

# Plot Z2 comparison
plot!(p1, t_vals, z2_num, label="Numerical z2", lw=4, color=:lightpink)
plot!(p1, comparison_data.t, comparison_data.u[:, 2], 
      label="Analytical z2", lw=2, linestyle=:dash, color=:red)

display(p1)

# Verification check
err = maximum(abs.(z1_num .- comparison_data.u[:, 1]))
println("Maximum Error in Z1: $err")
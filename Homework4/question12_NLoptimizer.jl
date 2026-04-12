# ===== Engineering Mathematics -- Homework 4 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 4th April, 2026

using JuMP
using Ipopt

# -- Define Surfaces -- #
# Define optimization surface
z(x,y) = (x * (x-y))^2 - exp(-x^2-y^2) - cos(x)

# Define constraint plane
p(x,y) = -3*(x+y)


function myFmincon()
    # Init optimizer model with Ipopt solver
    model = Model(Ipopt.Optimizer)
    
    # Quiet!
    set_silent(model)

    # Define variables
    @variable(model, -1.0 <= x <= 1.0)
    @variable(model, -1.0 <= y <= 1.0)
    @variable(model, z)

    # Initial guesses
    set_start_value(x, 0.2)
    set_start_value(y, 0.3)
    set_start_value(z, -1.5)

    # Define objective
    @objective(model, Min, z)

    # Define constraints
    @constraint(model, plane, x + y + z/3 == 0)
    
    # The unit circle inequality constraint
    @constraint(model, circle, x^2 + y^2 <= 1.0)
    
    # The surface equality constraint
    @constraint(model, surface, z == (x * (x - y))^2 - exp(-x^2 - y^2) - cos(x))

    # Run!
    optimize!(model)

    return value(x), value(y), value(z), termination_status(model)
end

# Execute the function
x_opt, y_opt, z_opt, status = myFmincon()

println("Solver Status: ", status)
println("Optimal x ≈ ", round(x_opt, digits=3))
println("Optimal y ≈ ", round(y_opt, digits=3))
println("Optimal z ≈ ", round(z_opt, digits=3))
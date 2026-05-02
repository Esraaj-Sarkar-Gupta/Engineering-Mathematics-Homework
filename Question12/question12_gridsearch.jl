# ===== Engineering Mathematics -- Homework 4 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 4th April, 2026

# Here a constrained grid search is implemented.

# -- Define Surfaces -- #
# Define optimization surface
z(x,y) = (x * (x-y))^2 - exp(-x^2-y^2) - cos(x)

# Define constraint plane
p(x,y) = -3*(x+y)

function run_grid_search()
    step = 0.0001
    tol  = 0.001

    X_best = (0.0, 0.0)
    z_min = Inf64

    # Range limits
    r = 1.0
    x_range = -r:step:r
    y_range = -r:step:r

    # ---- Begin Grid Search ---- #
    for x in x_range
        for y in y_range

            if (x^2 + y^2) > r^2
                continue
            end

            z_surface = z(x,y)
            z_plane = p(x,y)

            if abs(z_plane - z_surface) < tol
                if z_plane < z_min
                    z_min = z_plane
                    X_best = (x,y)
                end
            end
        end
    end
    
    return X_best, z_min
end

@time X_best, z_min = run_grid_search()

println("Grid Search Results:")
println("x ≈ ", round(X_best[1], digits=3))
println("y ≈ ", round(X_best[2], digits=3))
println("z ≈ ", round(z_min, digits=3))




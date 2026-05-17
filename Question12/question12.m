% ===== Engineering Mathematics -- Homework 2, Question 7 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 20th March, 2026
% Introduction to Non-Linear Optimization
% ---- Grid Search ---- %

N = 2000;   % Grid Resolution (Reduced for stable memory allocation)

% X and Y values within a unit circle
x_grid = linspace(-1, 1, N);
y_grid = linspace(-1, 1, N);

% Meshgrid (Fixed 'y' to 'Y')
[X, Y] = meshgrid(x_grid, y_grid);

tic
% Compute the surface and plane geometries (Standardized to uppercase Z)
Z_surface = (X.*(X - Y)).^2 - exp(-X.^2 - Y.^2) - cos(X);
Z_plane = -3 * (X + Y);

% -- Apply Constraints -- %
inside_circle = (X.^2 + Y.^2 <= 1);
intersection_tolerance = 0.05;

% Fixed casing: Z_surf -> Z_surface, Z_plane matched
on_intersection = abs(Z_surface - Z_plane) < intersection_tolerance;

% "valid points" mask (Fixed typo: inside_cirlce -> inside_circle)
valid_points = inside_circle & on_intersection;

% Extract the minimum Z from only the valid points -- apply the mask
Z_valid = Z_plane;
Z_valid(~valid_points) = NaN; % Turn invalid points to NaN so min() ignores them


[z_min_grid, min_idx] = min(Z_valid(:));
toc

% Handle edge case where no points satisfy the constraints
if isnan(z_min_grid)
    error('No valid points found within the intersection tolerance. Try increasing intersection_tolerance.');
else
    [row, col] = ind2sub(size(Z_valid), min_idx);
    x_grid_opt = X(row, col);
    y_grid_opt = Y(row, col);
    
    fprintf('Grid Search Minimum: z = %.4f at (x,y) = (%.4f, %.4f)\n', z_min_grid, x_grid_opt, y_grid_opt);
end
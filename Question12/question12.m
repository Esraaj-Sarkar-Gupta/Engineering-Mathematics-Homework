% ===== Engineering Mathematics -- Homework 2, Question 7 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 20th March, 2026

% Introduction to Non-Linear Optimization

% ---- Grid Search ---- %
N = 10000   % Grid Resolution

% X and Y values within a unit circle
x_grid = linspace(-1, 1, N);
y_grid = linspace(-1, 1, N);
;
% Meshgrid
[X, y] = meshgrid(x_grid, y_grid);

% Compute the surface and plane geometries
z_surface = (X.*(X - Y)).^2 - exp(-X.^2 - Y.^2) - cos(X);
z_plane = -3 * (X + Y);

% -- Apply Constraints -- %
inside_circle = (X.^2 + Y.^2 <= 1);

intersection_tolerance = 0.05;
on_intersection = abs(Z_surf - Z_plane) < intersection_tolerance;

% "valid points" mask
valid_points = inside_cirlce & on_intersection

% Extract the minimum Z from only the valid points -- apply the mask
Z_valid = Z_plane;
Z_valid(~valid_points) = NaN; % Turn invalid points to NaN so min() ignores them

[z_min_grid, min_idx] = min(Z_valid(:));
[row, col] = ind2sub(size(Z_valid), min_idx);
x_grid_opt = X(row, col);
y_grid_opt = Y(row, col);

fprintf('Grid Search Minimum: z = %.4f at (x,y) = (%.4f, %.4f)\n', z_min_grid, x_grid_opt, y_grid_opt);
% ===== Engineering Mathematics -- Homework 1 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 1st February, 2026

% -- Some commonly used helper functions -- %
function err = Error(x,y)
arguments
    x (:,1) double
    y (:,1) double = 0
end
    err = abs(x-y);
end

% ---- Question 2 ---- %

function dsin = dsin_numerical(x, epsilon)
arguments
    x (1,1) double = 1          % Accept coloumn vectors
    epsilon (1 ,1) double = 1e-3 % Default epsilon value
end

    dsin = (sin(x + epsilon) - sin(x)) ./ epsilon;
end

% -- 2 Part A: Compute dsinx(x = 1) -- The numerical derivative -- %

x_qn2 = 1.0;
numerical_derivative_sinx = dsin_numerical(1.0);
analytical_derivative_sinx = cos(x_qn2);

fprintf("The numerically computed derivative is %.5f and the analytically computed derivative is %.5f\n", ...
    numerical_derivative_sinx, analytical_derivative_sinx);

error = Error(numerical_derivative_sinx, analytical_derivative_sinx);

fprintf("The error is %.5f\n", error)

% -- 2 Part B : Generating the h vs err graph -- %

epsilon_values = logspace(-1, -16, 100);
numerical_derivative_values = arrayfun(@(e) dsin_numerical(x_qn2, e), epsilon_values);

error_arr = Error(numerical_derivative_values, analytical_derivative_sinx);

% Plot graph
loglog(epsilon_values, error_arr);
grid on
xlabel('\epsilon');
ylabel('|error|');
title("\epsilon vs Error graph");
shg

%{
    From the graph it is evident that the best value for epsilon lies
    around the value 1e-8, after which roundoff errors pick up.
%}

% -- 2 Part C : Using eps to estimate ideal h (hoipefully for extra credit) -- %

h_eps = sqrt(eps(max(1, abs(x_qn2)))) * max(1, abs(x_qn2)); % Use matlab eps

d_eps = dsin_numerical(x_qn2, h_eps);
err_eps = Error(d_eps, cos(x_qn2));

fprintf("h from eps is %.5e\n", h_eps);
fprintf("Derivative using h_eps is %.5f, error is %.5e\n", d_eps, err_eps);

% Compare with sweep minimum
[err_min, idx] = min(error_arr);
h_min = epsilon_values(idx);

fprintf("Best h from sweep is %.5e with error %.5e\n", h_min, err_min);

hold on
loglog(h_eps, err_eps, 'o')
loglog(h_min, err_min, 'x')
legend('|error|', 'h\_eps', 'h\_min', 'Location', 'best')
hold off

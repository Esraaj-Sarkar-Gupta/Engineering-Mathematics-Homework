% ===== Engineering Mathematics -- Homework 2, Question 6 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 8th February, 2026

% Housekeeping
clc; clear; close all;

% ==== Higher Order Methods ==== %

%% Part A -- Differentiation using Central Difference %%


% -- Function to be Differentiated -- %
f = @(x) sin(x);

% True analytic derivative
df = @(x) cos(x);

% ---- Numerical Evaluation -- Central Difference Method ---- %
%{
    delta_f(x) ~ f(x + h) - f(x-h) / 2h
%}

X0  = 1.0;
%{
    The point 0 is particularly bad for this analysis since
    sin(h) - sin(-h) ~ 2h.

    Since the program returns a 0, and thus error = 0 - 0, the plot
    attempts to plot log(0 - 0 = 0), which is undefined and thus it cannot
    plot.
%}

% -- List of h values -- %
h_values = logspace(-2, -16, 100);
evaluation = central_difference(f, X0, h_values);

% ---- Error Evaluation ---- %
true_val = df(X0);
errors = abs(true_val * ones(size(h_values)) - evaluation);

% Plot
figure;
loglog(h_values, errors, 'o-', 'LineWidth', 1.6, 'MarkerSize', 4);
xlabel('h', 'FontSize', 12)
ylabel('Absolute error', 'FontSize', 12)
title('Central difference error for d/dx(sin x) at x=0', 'FontSize', 13)
legend('Central difference', 'Location', 'best')
set(gca, 'FontSize', 11, 'LineWidth', 1.1)
hold off
shg;

% ---- Comparing h error vs Forward Difference Method ---- %
evaluation_forward_difference = forward_difference(f, X0, h_values);
errors_forward_difference = abs(true_val * ones(size(h_values)) - evaluation_forward_difference);

% Comparision plot
figure;
loglog(h_values, errors, 'o-', 'LineWidth', 1.6, 'MarkerSize', 4);
hold on
loglog(h_values, errors_forward_difference, 's--', 'LineWidth', 1.6, 'MarkerSize', 4);
hold off
xlabel('h', 'FontSize', 12)
ylabel('Absolute error', 'FontSize', 12)
title('Comparison between central difference and forward difference methods for sin(x = 0).', 'FontSize', 13)
legend('Central difference', 'Forward difference', 'Location', 'best')
set(gca, 'FontSize', 11, 'LineWidth', 1.1)
shg;


% ---- Functions ----- #
function delta_f = central_difference(f, x, h)
    delta_f = ( f(x + h) - f(x - h)) ./ (2*h);
end

function delta_f = forward_difference(f, x, h)
    delta_f = ( f(x + h) - f(x) ) ./ h;
end


%{
    OBSERVATIONS:
        > Optimal h for central diff is much larger (~ x100) than forward
        diff. However the abs error for central diff is much lower.

        > The linear parts of the graph make it clear that the slope of the
        central diff curve is greater than the forward diff graph. This
        observation agrees with what one would expect from the mathematical
        expressions of both these numerical methods.

        > We see that the error of the central difference method decreases
        with h^2. We may also show this to be true using the Taylor
        expansion of the series.
%}


%% Part B -- Solving ODEs using the Midpoint Method %%

% Housekeeping
clear;

% -- System -- %
p = [1];     % Parameters
X0 = [4];    % Intial conditions
h = 1e-4;    % Step size
T = 5;       % Final time

% -- Main Function -- %
[t, X] =  RK2(X0, p, h, T);

% Plotting
figure;
plot(t, X)
xlabel('t'); ylabel('X(t)')
title('RK2 (Midpoint Method) solution')
grid on
shg

% -- Local Functions -- %
function dX = g(X, t, p)
arguments
    X
    t
    p (1,1) double = 1
end
    % Example: time-dependent system dX/dt = p * t * X
    % (uses t explicitly; p is correctly applied as a scalar parameter)
    dX = p * t .* X;
end


function [t, X] = RK2(X0, p, h, T)
    N = round(T/(2*h)) + 1;     % Number of steps
    t = linspace(0, T, N).';    % Time values
    X = zeros(N, 1);            % Init X values
    X(1) = X0;                  % Initial value

    for k = 1:N-1
        X_mid = X(k) + g(X(k), t(k), p) * (h/2);
        X(k+1) =  X(k) + h * g(X_mid, t(k) + h/2);
    end
end
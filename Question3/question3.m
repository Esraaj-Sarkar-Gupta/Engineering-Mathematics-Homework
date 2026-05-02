% ===== Engineering Mathematics -- Homework 1 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 1st February, 2026

% Housekeeping
clear; clc; close all;

% ---- Question 3 -- Euler: error vs step size ---- %
% ODE: dx(t) = p x(t),  x(0)=X0

p  = -1;     % Parameter
X0 =  4;     % Initial condition
T  =  5;     % Final time

% Step sizes
h_values = logspace(-4, -7, 100).';

% Exact value at T (solved analytically by esraaj)
xT_exact = X0 * exp(p*T);

% Arrays to store resultant data
err_true_T   = zeros(size(h_values));   % |x_h(T) - x_true(T)|

% Relative error esitmates
err_est_T    = zeros(size(h_values));   % |x_h(T) - x_h/2(T)|
err_est_2x   = zeros(size(h_values));   % 2*|x_h(T) - x_h/2(T)|

% --- Main loop-- sweep h --- %
for i = 1:numel(h_values)
    h = h_values(i);

    % Euler with step h
    [t_h, x_h] = Eulers_Method(X0, p, h, T);
    xT_h = x_h(end);

    [t_h2, x_h2] = Eulers_Method(X0, p, h/2, T);
    xT_h2 = x_h2(end);

    err_true_T(i) = abs(xT_h - xT_exact);

    err_est_T(i)  = abs(xT_h - xT_h2);

    err_est_2x(i) = 2 * err_est_T(i);
end

% ---- Plot: true error and the estimate ---- %
figure('Name','Euler: error vs h','NumberTitle','off');
loglog(h_values, err_true_T, 'LineWidth', 2); hold on;
loglog(h_values, err_est_T,  'LineWidth', 2);
loglog(h_values, err_est_2x, 'LineWidth', 2);
grid on
xlabel('h');
ylabel('Error at t = T');
title('Euler method: true error vs step-doubling estimate');
legend('|x_h(T) - x_{exact}(T)|', '|x_h(T) - x_{h/2}(T)|', '2|x_h(T) - x_{h/2}(T)|', ...
       'Location','best');
shg

%{
Observations:
One may estimate the true error of a system using this method of comparing
the difference between the results of consecutive step sizes.

This is possible since error is directly propotionate to the step size --
euler is linear.

x_h(T) - x_h/2 = Ch - Ch/2 + O(h^2) + ... (we ignore non-linear terms).
Thus, the estimated error is offset by some constant from the true error by
a linear approximation.
%}

% -- Local Functions -- %

function dX = f(X, p)
    dX = p * X;
end

function [t, X] = Eulers_Method(X0, p, h, T)
    % Euler method for scalar ODE X' = f(X,p), with X(0) = X0
    % We force the grid to land exactly on T.

    N = round(T/h) + 1;       % number of points
    t = linspace(0, T, N).';  % time vector
    h = t(2) - t(1);          % actual step used

    X = zeros(N,1);
    X(1) = X0;

    for k = 1:N-1
        X(k+1) = X(k) + h * f(X(k), p);
    end
end

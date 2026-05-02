% ===== Engineering Mathematics -- Homework 1 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 1st February, 2026

% ---- Question 4 -- Euler's Method ---- %

p = [-1];   % Parameter list
X0 = [4];   % Initial Conditions: x_0,
h = 0.01;   % Step size
T = 5;      % Final time

% Main Function
[t, X] = Eulers_Method(X0, p, h, T);

% Plotting
plot(t, X)
xlabel('t'); ylabel('X(t)')
title('Euler method solution')
grid on
shg

% -- Local Functions -- %

function dX = f(X, p)
    dX = p * X;
end

function [t, X] = Eulers_Method(X0, p, h, T)
    N = round(T/h) + 1;         % Number of steps
    t = linspace(0, T, N).';    % Time values
    X = zeros(N,1);             % X values
    X(1) = X0;                  % Initial X value
    
    % Main loop
    for k = 1:N-1
        X(k+1) = X(k) + h * f(X(k), p);
    end
end
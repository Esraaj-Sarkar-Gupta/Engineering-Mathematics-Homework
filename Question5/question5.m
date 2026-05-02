% ===== Engineering Mathematics -- Homework 1 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 1st February, 2026

% Housekeeping
clear; clc; close all;

% ---- Question 5 -- Equations of Motion of a Double Pendulum ---- %

% -- System Parameters -- %
g = 9.8;
m1 = 2;
m2 = 3;
R = 10;
delta = 3;

% Simplify parameteras (for the sake of our sanity)
A = (m1 * R^2 / 3) + (1.5 * m2 * delta^2);
B = 1.5 * delta^2;
C = 2 * m2 * g * R;
D = m2 * B;
E = 2 * m2 * g * delta;

% Parameter Vector
p = [A; B; C; D; E];

% -- Time Spans -- %
t0 = 0;
T  = 20;

% Euler step size
h = 1e-6;

% -- Intial Conditions -- %
y0 = [1.0; 0.0; 0.0; 0.0]; % y = [theta, theta_dot, phi, phi_dot]

% -- Main Function --%
rhs = @(t,y) coupled_system_rhs(t, y, p);

% Solve using Euler
[t, Y] = euler_system(rhs, t0, T, y0, h);

theta     = Y(:,1);
theta_dot = Y(:,2);
phi       = Y(:,3);
phi_dot   = Y(:,4);

% ---- Plotting ---- %
figure('Name','Angles vs time','NumberTitle','off');
plot(t, theta, 'LineWidth', 2); hold on;
plot(t, phi,   'LineWidth', 2);
xlabel('Time');
ylabel('Angle (radians)');
legend('\theta_1(t)','\theta_2(t)','Location','best');
grid on;
title('Euler Method Solution of the Coupled System');
shg

figure('Name','Phase space','NumberTitle','off');

subplot(1,2,1);
plot(theta, theta_dot, 'LineWidth', 1.5);
xlabel('\theta_1'); ylabel('\dot{\theta_1}');
title('Phase Space: \theta_1 vs. \dot{\theta_1}');
grid on;

subplot(1,2,2);
plot(phi, phi_dot, 'LineWidth', 1.5);
xlabel('\theta_2'); ylabel('\dot{\theta_2}');
title('Phase Space: \theta_2 vs. \dot{\theta_2}');
grid on;

shg

%{
    One may observe non-linear patterns and chaotic aperiodicity in these
    graphs, indicative of a double pendulum's chaotic motion. This is
    especially evident in the non-closed phase diagrams.
%}

% == Local Functions == %
%{
    Likely to Professor Ruina's disappointment, I derived these equations
    of motion of a double pendulum using Lagrangian Mechanics as opposed to
    having used moment balance.

    One day I will learn to do better.
%}

function dydt = coupled_system_rhs(~, y, p)
    % Unpack parameters
    A = p(1); B = p(2); C = p(3); D = p(4); E = p(5);

    % y = [theta, theta_dot, phi, phi_dot]
    theta     = y(1);
    theta_dot = y(2);
    phi       = y(3);
    phi_dot   = y(4);

    % Equations of Motion
    theta_ddot = (B * E * cos(phi) - C * D * sin(theta)) / (A * D);
    phi_ddot   = (-A * E * cos(phi) - B * E * cos(phi) + C * D * sin(theta)) / (A * D);

    dydt = [theta_dot; theta_ddot; phi_dot; phi_ddot];
end

function [t, Y] = euler_system(rhs, t0, tf, y0, h)
    % Euler integrator for vector ODEs y' = rhs(t,y)
    arguments
        rhs (1,1) function_handle
        t0  (1,1) double
        tf  (1,1) double
        y0  (:,1) double
        h   (1,1) double {mustBePositive}
    end

    N = floor((tf - t0)/h) + 1;
    t = t0 + (0:N-1)'*h;

    nstate = numel(y0);
    Y = zeros(N, nstate);
    Y(1,:) = y0.';

    for k = 1:N-1
        yk = Y(k,:).';
        Y(k+1,:) = (yk + h * rhs(t(k), yk)).';
        % If you start seeing NaNs, your step is too big (Euler is fragile).
    end
end

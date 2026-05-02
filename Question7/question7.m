% ===== Engineering Mathematics -- Homework 2, Question 7 ===== %
% Author: Esraaj Sarkar Gupta
% Date: 8th February, 2026

% Housekeeping
clc; clear; close all;

%% Part A -- Equations of Motion %%

%{
    Equations of Motion (derived analytically)

    x_ddot = (1/m) [ -k ( sqrt(x^2 + y^2) - L0) \hat{e} dot i]
    y_ddot = (1/m) [(-mg) - k ( sqrt(x^2 + y^2) - L0 \hat{e} dot j]
    
        where \hat{e} = (xi + yj) / sqrt(x^2 + y^2)
%}

%% Part B -- Plotting the Trajectory of the Mass

% System parameter struct
params.m  = 1.0;
params.k  = 1.0;
params.L0 = 0.0;
params.g  = 9.81* 0;

S0 = [1.0; 0.0; 0.0; 1.0];          % [x0; y0; vx0; vy0]
tspan = [0 10];                      % Initial time, Final time

h = 1e-3;

% -- SOLVE -- %
[t, S] = RK2(@(t,S) system(t,S,params), tspan, S0, h);

% -- Plotting Results -- %

x = S(:,1);
y = S(:,2);

figure;
plot(x, y, 'LineWidth', 1.5);
axis equal
grid on
xlabel('x'); ylabel('y');
title('Trajectory of the mass (RK2 Midpoint)');

% -- Animation -- %
fig = figure;
ax = axes(fig);
cla(ax, 'reset');
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x'); ylabel(ax, 'y');
title(ax, 'Animation of System');

% -- Lock Figure Axes -- %
pad = 0.1;  % 10% padding
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);

xr = xmax - xmin;  yr = ymax - ymin;

if xr == 0, xr = 10; end
if yr == 0, yr = 10; end

xlim(ax, [xmin - pad*xr, xmax + pad*xr]);
ylim(ax, [ymin - pad*yr, ymax + pad*yr]);

axis(ax, 'manual');

% Create handles
mass_pt   = plot(ax, x(1), y(1), 'o', 'MarkerSize', 10, 'LineWidth', 2);
spring_ln = plot(ax, [0 x(1)], [0 y(1)], '-', 'LineWidth', 1.5);
trail_ln  = plot(ax, x(1), y(1), '-', 'LineWidth', 1.2);

trail_len = 300;
skip = 5;

for k = 1:skip:length(t)
    % -- Dynamically update plot elements -- %
    set(mass_pt,   'XData', x(k),     'YData', y(k));
    set(spring_ln, 'XData', [0 x(k)], 'YData', [0 y(k)]);

    k0 = max(1, k - trail_len);
    set(trail_ln,  'XData', x(k0:k),  'YData', y(k0:k));

    drawnow;
end


% -- Local Functions -- %

% Modelling the System
function dS = system(t, S, params)

% -- Unpack State -- %
% Python tuple unpacking is better <\3
x  = S(1);
y  = S(2);
vx = S(3);
vy = S(4);

% -- Unpack parameter structure -- % 
m  = params.m;
k  = params.k;
L0 = params.L0;
g  = params.g;

% -- Geometry -- %
r2 = x^2 + y^2;
r  = sqrt(r2);
if r == 0
    e_x = 0;
    e_y = 0;
else
    e_x = x / r;
    e_y = y / r;
end

deltaL = r - L0;

% -- Forces -- %

F_spring_mag = -k * deltaL;
F_spring_x   = F_spring_mag * e_x;
F_spring_y   = F_spring_mag * e_y;

F_grav_x = 0;
F_grav_y = -m * g;

F_net_x = F_spring_x + F_grav_x;
F_net_y = F_spring_y + F_grav_y;

ax = F_net_x / m;
ay = F_net_y / m;

% -- Pack Derivatives -- %

dx  = vx;
dy  = vy;
dvx = ax;
dvy = ay;


dS = [dx; dy; dvx; dvy];
end

function [t, S] = RK2(f, tspan, S0, h)
    t0 = tspan(1);
    tf = tspan(2);

    N = floor((tf - t0)/h) + 1;      % Number of steps
    t = (t0 + (0:N-1)*h).';          % Time values

    S = zeros(N, numel(S0));         % State history
    S(1,:) = S0(:).';                % Initial state

    for k = 1:N-1
        tk = t(k);
        Sk = S(k,:).';               % column state

        k1 = f(tk, Sk);
        S_mid = Sk + (h/2)*k1;

        k2 = f(tk + h/2, S_mid);
        S_next = Sk + h*k2;

        S(k+1,:) = S_next.';
    end
end

%% Part C -- Checking the Numerical Solution

%{ 
    i. No spring (k=0). One would assume that the mass simply falls, given
    that the initial velocities are 0. With non-zero initial velocities one
    must expect to observe parabolic flight.
        > We verify these cases to be true.
    
    ii. Rigid Connection (k = 10e10). This will make the spring behave like
    a rigid connection. This should ccause the system to behave like a
    pendulum.
        > Large values result in numerical errors. Marginally large values
        (approx 1e3) leads to unintuitive behavior.
    
    iii. g = 0 and \vec{v} = 0 must constrain the motion of the mass to a
    straight line.
        > From the animation we see this to be the case.
    
    iv. When L0 = 0, in the equations of motion, we see that the spring
    force loses its mag(r) contribution altogether. This makes the force
    perfectly linear along both axes, and thus the motion of the mass must
    follow a radial line. This is the equation of a harmonic osciallator
    in a 2D space.
        > From the numerical results we may verify this to be true.

%}

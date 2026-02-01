%Euler's method
% Version  2    
% Andy Ruina  29 January 2026
% NEW THINGS
%   * coupled ODEs  (P & S). State z = [P,S]'
%   * RHS in a separate matlab function: myRHS
%   * ODE solver in a separate function: ODEEuler
%   * matlab struct p for parameters
%   * passing a function name ('handle') to another function using @

% Much of the structure of the data in this example is chosen
% to match the syntax used by Matlab's ODE solvers.
% For example, one can't directly pass parameters into the solver.
% This is avoided using 'anonymous' functions and the @ symbol

%HEADER STUFF
% A common error is that old work is contaminating your new work.
% Or, an old plot shows and you think it is new.  To avoid this, use 
% the commands below at the top of script files.
clear; close all % clears all variables and figures. Good practice.

disp('-----------------------------------')
disp(['Start at '  datestr(now)])
tic     % Start the stop watch. Used with toc to time the program.


%Number of points in the solution array
n = 1e3 +1;

% Initial and final time
t0 = 0; tfin = 6*pi;
tspan = [t0 tfin];   % pack t0 and tfin into one vector

% Initial conditions 
P0  = 1;  S0  = 0;
z0  = [P0;S0];        % z is the state, consisting of P and S

% Set ODE parmeters
p.cP = 1;
p.cS = 1;


% SOLVE ODES
%
% First, create a new RHS function that does not have parameters.
% The parameters in myRHS are frozen at their values at the time
% of the execution of this line of code.
% Using @ for a function handle is a bit confusing. Read about 
% it on GPT or whatever. The next line creates a new function RHS 
% from myRHS, that is a function only of t and z with the parameters 
% p frozen.

RHS = @(t,z)  myRHS(t,z,p);  % gets rid of the p which the solver cannot accept.

% Here is where we call the ODE solver, in this case ODEEuler.
% It is below.  Later we will call professional solvers.
% It needs the name of the RHS file, the time interval and the ICs
% The RHS function whose 'handle' (name) is passed can only have 
% the arguments t and z.  It can't have p.  Thus the use of the 
% anonymous function above.

    [tarray zarray] = ODEEuler(RHS, tspan, z0, n);   %Command to solve ODEs

% At this point the ODE has been solved.  The rest is just display and
% plotting


% Unpack ODE solution
   Parray = zarray(:,1);  % P vs time is first   col of zarray
   Sarray = zarray(:,2);  % S vs time is second  col of zarray


%  PLOTING
%%%%%%
figure(1)
plot(tarray, Parray, tarray, Sarray, LineWidth= 4)
title('Phil and Sally feelings vs time')
xlabel('time');  ylabel('P and S''s feelings');
legend('Phil', 'Sally')
shg   %show graph

%%%%%
figure(2)
plot(Parray, Sarray)
title('Trajectory in state space')
xlabel('Phil''s feelings')
ylabel('Sally''s feelings')
axis equal
shg


% The stuff below helps keep track of things.
disp(['Done at  '  datestr(now)])
toc 
disp(['Number of time steps = ' num2str(n)])
disp('-----------------------------------')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%    myRHS
%    ODEEuler

% DEFINE RHS  (This could be in a separate .m file).
%%%%%%%%%%%%%%%




z
 for k = 1:n-1
   zold = znew;    % the former new state is now the old state
   t = tarray(k);  % present time

   % Euler's method:  the core, the heart, of this code is
   % the set of two lines below. Everything else is built around that.
   zdot = funct(t,zold);       %  This is the ODE
   znew = zold + zdot * h;     %  One step of Euler's method

   zarray(k+1,:) = znew';      %  Save the state as a row of zarray    
  end  % End of Euler method loop
 end   % End of Euler Method funcion
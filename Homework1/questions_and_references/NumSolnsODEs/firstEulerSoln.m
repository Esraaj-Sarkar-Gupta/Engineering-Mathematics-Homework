%Euler's method
% Version  1     
% Andy Ruina   - 27 January 2026

% Initial and final time
t0 = 0; tfin = 5;
z0 =  4;
%Number of points in the solution array
n = 1001;
%initialize tarray and zarray
tarray =  linspace(t0, tfin, n);
zarray = zeros(n,1);
zarray(1)   = z0;

h = (tfin - t0)/n;

%RHS function
tic
for k = 1:n-1
  zold = zarray(k);
  t = tarray(k);

  %Euler's method
  f = -zold;
  znew = zold + f* h;

  zarray(k+1) = znew;
    
end
toc

plot(tarray, zarray)
title('ODe solution')
xlabel('time');  ylabel('z');
shg   %show graph

disp('Calculation is done.')


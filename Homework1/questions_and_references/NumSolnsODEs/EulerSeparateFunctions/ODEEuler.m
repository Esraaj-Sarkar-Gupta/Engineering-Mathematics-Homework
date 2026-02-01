 function  [tarray, zarray] = ODEEuler(funct, timespan, z0, n)
 % General purpose ODE solver for systems of 1st order equations
 %    funct    =  handle (name of) of RHS file
 %    timespan =  [t0, tfin]
 %    z0       =  ICs
 %    n        =  Number of time points = # steps + 1

 %Unpack inputs
 t0    = timespan(1); tfin = timespan(2); % unpack tspan
 neqns = length(z0);       % Number of state vars = Num of eqns

 h     = (tfin - t0)/(n-1);% h= Delta t = duration of time step
   
 % Initialize tarray and zarray
 tarray = linspace(t0,tfin, n);
 zarray = zeros(n,neqns);  %  one row for each time step
 zarray(1,:)  = z0';       %  IC put in first row  of output
 znew         = z0 ;       %  Used for the first execution of loop below


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
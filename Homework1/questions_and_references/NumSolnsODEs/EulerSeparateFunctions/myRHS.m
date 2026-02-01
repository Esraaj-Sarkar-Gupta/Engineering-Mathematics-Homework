function  zdot = myRHS(t,z,p)  % This has parameters p
% Given the state z and time t,
% calculate zdot (the rate of change of z)
% This RHS doesn't happen to use t, but we include it 
% because the solver passes it.

% Unpack t, z, and p
  cP = p.cP;  cS = p.cS;     % unpack the struct p

% Unpack z
  P = z(1);     % Unpack z
  S = z(2);

% The ODEs, The rules governing evolution of Phil and Sally's emotions
  Pdot =  - cP * S ;
  Sdot =  + cS * P;

% Pack readable variables back into zdot
    zdot = [Pdot; Sdot];   %column vector is needed
end % end of myRHS function

% Note:  This whole function (myRHS) could be written with one line of code:
%         zdot = [ -z(1) * p.cP;  z(2) * p.cS];
% But, this is not as  readable as the Pdot & Sdot eqs above.
% This issue is exacerbated for more complex problems.
% I recommend writing ODEs with readable variables (like P and S) and
% not with elements of z (like z(1) and z(2)).
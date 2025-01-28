
function [uSS,P] = SSprofile(p,P0,RelTol)
%% [uSS,P] = SSprofile(p,P0,RelTol)
% Returns the steady-state profile of the states for a given downhole
% pressure
% 
% Syntax 
% 
% [uSS,P] = SSprofile(p,P0) computes the steady-state profiles of the
% conservative states (n,m,I), vertically concatenated in the uSS vector,
% for a given value of the bottom hole Pressure P0 and a parameter struct
% p.
% 
% [uSS,P] = SSprofile(p,P0,RelTol) does the same, with an accuracy defined
% by the RelTol parameter. Default value for RelTol is 1e-6.
% 
% Description
% 
% This function computes the mass inflow rates of gas and liquid from the
% bottom pressure, and solves the steady-state equation for pressure
% along the well. The conservative states are then computed using the
% FluxPressure2state.m function. When the bottom pressure P0 provided by
% the user is too low and does not provide a viable solution, an error is
% returned. 
%
% 01.04.2014    Ulf Jakob F. Aarsnes & F. Di Meglio


if nargin < 3
    RelTol = 1e-6;
end

[W_Gres,W_Lres,~,~] = LeftMassrates(p,P0);

I_G0 = (W_Gres + p.W_Gbit)/p.A;
I_L0 = (W_Lres + p.Q_Lbit.*p.rho0_L)/p.A;

options = odeset('Events',@(s,P,I_G,I_L,p)events(s,P,I_G,I_L,p),'RelTol',RelTol);
[~,P,TE]=...
    ode45(@odeSSfun,p.s,P0,options,I_G0,I_L0,p);

if  ~isempty(TE)
    % Integration did not converge.
    error('cmpSS:PTendsTo0','Integration did not converge and was terminated.')
end
I_G = I_G0*ones(p.P,1);
I_L = I_L0*ones(p.P,1);
[n,m,I] = FluxPressure2States(p,I_G,I_L,P);
uSS = u2u(p,n,m,I);
end


function [value,isterminal,direction] = events(s,pressure,I_G,I_L,p)
% This function ensures that the solver of the steadt-state equation does
% not try to solve equations for which there are no solutions. More
% precisely, when the computed pressure is about to become negative, the
% solver should stop integration and return an error. 

Pdot = odeSSfun(s,pressure,I_G,I_L,p);

if s==0
    value = 1;
else
    value = pressure + .8*Pdot*(p.L-s);
end
isterminal = 1;   % stop the integration
direction = 0;   % Direction doesn't matter

end





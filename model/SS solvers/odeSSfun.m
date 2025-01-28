function dPds = odeSSfun(s,pressure,I_G,I_L,p)
%% function dPds = odeSSfun(s,pressure,I_G,I_L,p)
%
% Returns the steady-state pressure gradient at position 's'
% 
% Syntax
% 
% dPds = odeSSfun(s,pressure,I_G,I_L,p) returns the pressure gradient at
% the spatial position s, as a function of the pressure at s, and the mass
% influxes I_G and I_L. 
% 
% Description
% 
% This function basically computes the right-hand-side of the steady-state
% pressure equations P'(s) = f(P(s),I_G,I_L). It is meant to be fed to an
% ODE solver like ode45. 
% c_G=(interp1q(p.s,p.c_G,s));

c_G = p.c_G(1)+(s-p.s(1))/(p.s(end)-p.s(1))*(p.c_G(end)-p.c_G(1));

[n,m,I] = FluxPressure2States(p,I_G,I_L,pressure,1,s);

[~, v_G, ~, ~, alpha_L, F_W, F_G, v_M,...
    rho_G,rho_L,~,Delta] =...
    variablesFromStates(p,n,m,I,1,s);


dvMdP = -I_G/c_G^2/rho_G^2-I_L/p.c_L^2/rho_L^2;
dDeltadP = 2*p.K*dvMdP*(v_M*p.K+p.S)+4*I_G/c_G^2/rho_G^2*(v_M*(p.K-1)+p.S)-...
    4*I_G/rho_G*(dvMdP*(p.K-1));
dvGdP = ((p.K*dvMdP+dDeltadP/2/sqrt(Delta)))/2;
dalphaGdP = -I_G*(rho_G*dvGdP+v_G/c_G^2)/v_G^2/rho_G^2;
dvLdP = -I_L*(-rho_L*dalphaGdP+alpha_L/p.c_L^2)/alpha_L^2/rho_L^2;

dPds=(F_W+F_G)/(I_G*dvGdP+I_L*dvLdP+1);
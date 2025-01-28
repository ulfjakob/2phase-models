function [n,m,I] = FluxPressure2States(p,I_G,I_L,pressure,k,s)
% FluxPressure2States
% 
% Computes the conservative states when given influxes and pressure
% distribution
% 
% SYNTAX
% 
% [n,m,I] = FluxPressure2States(p,I_G,I_L,pressure) computes the
% conservative states n, m and I over the whole spatial domain.
%  --Inputs--
% 
%  - p: parameter struct
%  - I_G (scalar): mass influx of gas [kg/s/m2]
%  - I_L (scalar): mass influx of liquid [kg/s/m2]
%  - pressure (vector): pressure distribution along the well [Pa]
%  
% [n,m,I] = FluxPressure2States(p,I_G,I_L,pressure,k) computes the same,
% but only for the spatial indices defined by k.
% 



if nargin == 4
    k = 1:p.P;
    c_G = p.c_G(k);
elseif nargin ==6
    c_G = p.c_G(1)+(s-p.s(1))/(p.s(end)-p.s(1))*(p.c_G(end)-p.c_G(1));
elseif nargin == 5;
    c_G = p.c_G(k);
end


rho_G = pressure./c_G.^2;
rho_L = pressure./p.c_L.^2+p.rho0_L;     %Compressible liquid

v_M = I_G./rho_G+I_L./rho_L;    %Mixture velocity

Delta = (v_M*p.K+p.S).^2 - 4*I_G./rho_G.*(v_M*(p.K-1)+p.S);

switch p.slipLaw
    case 'Evje'
        v_G = ((v_M*p.K+p.S)+sqrt(Delta))/2;% The other root is discarded (non-physical)
    otherwise
        error('SS solutions for other sliplaws than Evje have not been implemented');
end

alpha_G = I_G./v_G./rho_G;
alpha_L = 1-alpha_G;

n = alpha_G.*rho_G;
m = alpha_L.*rho_L;
I = I_G+I_L;

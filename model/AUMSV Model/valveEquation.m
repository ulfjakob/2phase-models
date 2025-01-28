function [W,gamma] = valveEquation(p,n,m,P,v_G,v_L,delta)
%% Function valveEquation
% 
% Syntax
% 
% W = valveEquation(p,n,m,P,v_G,v_L,delta) computes the mass flow
% rate W through a choke from the conservative variables n and m, the
% Upstream Choke Pressure P, the fluid velocities v_G and v_L and the
% regularizing parameter delta.
% 
% [W,gamma] = valveEquation(p,n,m,P,v_G,v_L,delta) also returns gamma,
% which is an approximation of the inverse of the square root of the
% mixture density. 
% 
% WARNING:  The choke model depends on the value of the parameter 
% p.valveModel, which can take the values 'Fjalstad', 'nonConstantY', or
% 'multiplier'.

I_G = n.*v_G;
I_L = m.*v_L;

% Liquid quality
chi_L=I_L/(I_G+I_L);

rho_G = P/p.c_G(end)^2;
rho_L = P/p.c_L^2+p.rho0_L;     %Compressible liquid


Phi_LO = 1;

switch p.valveModel
    case 'Fjalstad'
        Y = p.Y;
    case 'nonConstantY'
        X = min( (P-p.p_s)/P , p.x_PT );
        Y = (1-X*1.4/(3*p.gamma*p.x_PT))* X *(P/(P-p.p_s));
    case 'multiplier'
        % See [Schüller 2003]
        Phi_LO = sqrt(1 + (1-chi_L)*(rho_L/rho_G-1));
        Y = p.Y;
    otherwise
        W = nan;
        gamma = nan;
        return;
end

gamma = ( chi_L/sqrt(rho_L)+(1-chi_L)/sqrt(rho_G)/Y )*Phi_LO;

% Regularization of the square root. Type HELP regRoot for details
regSqrtP = regRoot(P-p.p_s,delta);

W = p.C_v*p.Z*regSqrtP/gamma;

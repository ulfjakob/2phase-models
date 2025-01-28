function [P, v_G, v_L, alpha_G, alpha_L, F_W, F_G, v_M,...
    rho_G, rho_L, Phi, Delta...
    dalpha_Gdq1,dalpha_Gdq2,dPdq1,dPdq2] =...
    variablesFromStates(p,n,m,I,k,s)
%% function [P,v_G,v_L,alpha_G,alpha_L,F_W,F_G,...
%   Phi,dalpha_Gdq1,dalpha_Gdq2,dPdq1,dPdq2] =...
%   variablesFromStates(p,n,m,I)
%
% Syntax
% 
% [P,v_G,v_L,alpha_G,alpha_L,F_W,F_G,...
%   Phi,dalpha_Gdq1,dalpha_Gdq2,dPdq1,dPdq2] =...
%   variablesFromStates(p,n,m,I)
% 
% computes the following physical variables:
% 
% - P: pressure [Pa]
% - v_G: gas velocity [m/s]
% - v_L: liquid velocity [m/s]
% - alpha_G: void fraction [m3/m3]
% - alpha_L: liquid volume fraction [m3/m3]
% - F_W: pressure loss due to wall friction (momentum equation) [Pa/m]
% - F_G: pressure loss due to gravity (momentum equation) [Pa/m]
% - Phi: Slip law function [m/s]
% - dalpha_Gdq1: partial derivative of alpha_G with respect to n
% - dalpha_Gdq2: partial derivative of alpha_G with respect to m
% - dPdq1: partial derivative of P with respect to n
% - dPdq2: partial derivative of P with respect to m
% 
% from the conservative variables (n,m,I), and using the parameters defined
% by the parameter struct p. 
% 
% Description
% This function is to be used whenever physical variables need to be
% computed from the conservative variables. 

if nargin == 5
    c_G = p.c_G(k);
    fVec = p.fVec(k);
    thetaVec = p.thetaVec(k);
elseif nargin ==6
    c_G = p.c_G(1)+(s-p.s(1))/(p.s(end)-p.s(1))*(p.c_G(end)-p.c_G(1));
    fVec = p.fVec(1)+(s-p.s(1))/(p.s(end)-p.s(1))*(p.fVec(end)-p.fVec(1));
    thetaVec = p.thetaVec(1)+(s-p.s(1))/...
    (p.s(end)-p.s(1))*(p.thetaVec(end)-p.thetaVec(1));

else
    c_G = p.c_G;
    fVec = p.fVec;
    thetaVec = p.thetaVec;
end

Mach = (c_G./p.c_L).^2;

% Solution of second order equation
delta = (p.rho0_L-Mach.*n-m).^2+4.*Mach.*n.*p.rho0_L;
[Delta,ddeltadDelta] = regMax0(delta,.01);

%Gas volume fraction
alpha_G = .5+(-Mach.*n-m+sqrt(Delta))./2./p.rho0_L;

alpha_L = 1-alpha_G;

P = (alpha_G<0.5) .* (m./max(alpha_L,.01)-p.rho0_L)*p.c_L^2 + ...
    (alpha_G>=0.5) .* n./max(alpha_G,.01).*c_G.^2;%The max( . ,0.1) are 
% only here to avoid the Nan. They do not impact the actual value of P.

% Liquid density
rho_L = p.rho0_L+P/p.c_L^2;

% Gas density
rho_G = P./c_G.^2;

% Slip law
switch p.slipLaw
    case 'Evje'
        %---Evje slip law---
        v_G = (I+p.S.*m./(p.K-(p.K-1).*alpha_G))./...
            (n+m.*(1-(p.K-1)./(p.K-(p.K-1).*alpha_G)));
        Phi = ((p.K-1).*v_G+p.S)./(p.K-(p.K-1).*alpha_G);
        v_L = v_G-Phi; 
    case 'ZF'
        % %---Zuber-Findlay slip law---
        if (alpha_L(1) > 1.d-10) % Mixture or Pure liquid
            v_G = (I + m.*p.S./(p.K.*(1-alpha_G)) ) ./...
                ( n + m.*(1 - (p.K-1)./(p.K.*(1-alpha_G))) );
            Phi = ((p.K-1)*v_G + p.S) ./ (p.K.*(1-alpha_G));
            v_L = v_G-Phi;
        else                  % Pure Gas
            v_G = (I + rho_L*p.S/p.K) ./ (n - rho_L.*(p.K-1)/p.K);
            v_L = 0.*v_G; % v_l non relevant
            Phi = 0.*v_G;
        end
    case 'Flatten'
        if (alpha_L(1) > 1.d-10) % Mixture or Pure liquid
            v_G = (I + m.*(.5*sqrt(1-alpha_G))./(p.K.*(1-alpha_G)) ) ./...
                ( n + m.*(1 - (p.K-1)./(p.K.*(1-alpha_G))) );
            Phi = ((p.K-1)*v_G + (.5*sqrt(1-alpha_G))) ./ (p.K.*(1-alpha_G));
            v_L = v_G-Phi;
        else                  % Pure Gas
            v_G = (I + rho_L.*(.5*sqrt(1-alpha_G))/p.K) ./ (n - rho_L.*(p.K-1)/p.K);
            v_L = 0.*v_G; % v_l non relevant
            Phi = 0.*v_G;
        end
end

% Mixture velocity
v_M = alpha_G.*v_G + alpha_L.*v_L;

% Friction source term
switch p.frictionModel
    case 'turbulent'
        F_W = -fVec.*(m+n).*abs(v_M).*v_M./(p.Dc-p.Dp);
    case 'Flatten'
        F_W = -32*v_M.*(alpha_G*p.eta_G+(1-alpha_G)*p.eta_L)./(p.Dc-p.Dp).^2;
    case 'linear'
        F_W = -fVec.*(m+n).*v_M./(p.Dc-p.Dp);
end

% Gravity source term
F_G = -(n+m).*p.g.*sin(thetaVec);


%% Jacobian

dDeltadq1 = -2*Mach.*(p.rho0_L-Mach.*n-m)+4.*Mach.*p.rho0_L;
dDeltadq2 = -2*(p.rho0_L-Mach.*n-m);

dalpha_Gdq1 = 1./2./p.rho0_L.*(-Mach+ddeltadDelta.*dDeltadq1./2./sqrt(Delta));
dalpha_Gdq2 = 1./2./p.rho0_L.*(-1+ddeltadDelta.*dDeltadq2./2./sqrt(Delta));

dPdq1 = (alpha_G>=0.5) .* (alpha_G-n.*dalpha_Gdq1)./max(alpha_G,.01).^2.*c_G.^2 ...
    + (alpha_G<0.5) .* m.*p.c_L.^2.*dalpha_Gdq1./max(alpha_L,.01).^2;

dPdq2 = (alpha_G>=0.5) .* (-n.*c_G.^2.*dalpha_Gdq2./max(alpha_G,.01).^2) ...
    + (alpha_G<0.5) .* (alpha_L+m.*dalpha_Gdq2)./max(alpha_L,.01).^2.*p.c_L.^2;










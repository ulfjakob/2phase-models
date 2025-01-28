% function [uu] = runAUMSVsteps(p, uu0, T0)
% Simulate the well model specified in p for p.dt seconds.
%
% Transient two-phase code based on AUSMV scheme: Gas and Water.
% Based on implementation by Steinar Evje.
%%
function [uu] = runAUMSVsteps(p, uu0, T0)
% Check if the struct is appropriate
if ~strcmp(p.valveModel,'constantP') && ~strcmp(p.valveModel,'Fjalstad') ... 
        && ~strcmp(p.valveModel,'MPD')
    error('Valve models has yet to be implemented for this simulator');
end

% Number of fluxes
Nf = p.P+1;
% Number of cells
N = p.P;
% Local timestep
% dt = 0.003;
% dt = p.ds/(1200);
dt = p.ds/2000;

dtdx = dt/p.ds;
t = T0;
Te = T0+p.dt; % Time for end of simulation
Nt = p.dt/dt;  %Number of total timesteps

% Initialize states
[n,m,I] = u2u(p,uu0);

% P1 = variablesFromStates(p,n(1),m(1),I(1),1);


% Check area: MERK HVORDAN VI KAN FORANDRE AREALET.
% Gemoetry below is 8.5" X 5" /typical 8 1/2" section well.
% do = 0.1*ones(p.P,1);
% di = 0*ones(p.P,1);
% areal = 3.14/4*(do.*do- di.*di);
% arear = 3.14/4*(do.*do- di.*di);
% area = 3.14/4*(do.*do- di.*di); 
% ang=3.14/2;
% The code below can be activated if one wants to introduce area changes inside
% the well.
%      for i = 25:p.P-1
%       do(i)=0.2;
%       di(i)=0.0;
%       areal(i+1) = 3.14/4*(do(i)*do(i)- di(i)*di(i));
%       arear(i) = 3.14/4*(do(i)*do(i)- di(i)*di(i));  
%       end 
%      arear(p.P) = arear(p.P-1);

% Initialize fluxes
flc = zeros(Nf,3);  % Liquid flux
fgc = zeros(Nf,3);  % Gas flux
fp = zeros(Nf,3);   % Pressure Flux

k = 1;
while t < Te
k = k+1;
t = t+dt; 
 
    % Trick to ensure numerical stability
    n = n.*(n>1e-8) + 1e-8.*(n<1e-8);
    m = m.*(m>1e-8);
      
    % Coefficients:
    a = 1/p.c_L.^2;
    b = p.rho0_L - m - p.c_G.^2.*n*a;
    c = -p.rho0_L*p.c_G.^2.*n;

    % Analytical solution:
    P = (-b+sqrt(b.^2 - 4*a.*c))/(2*a);  % Pressure
    rho_L = p.rho0_L + P/p.c_L^2;     % Density of liquid
    rho_G = P./p.c_G.^2;                   % Density of gas
    alpha_G = n./rho_G;
    alpha_L = 1-alpha_G;

    % Resolve slip law to compute velocities
    xint = (alpha_G-0.75)/0.25;
    k0 = p.K + (1.0*xint - p.K*xint).*(alpha_G >= 0.75).*(alpha_G <= 1);
    s0 = p.S + (0.0*xint - p.S*xint).*(alpha_G >= 0.75).*(alpha_G <= 1);
    
    k1 = 1*(alpha_G>=1-1e-6) + ...
        (1-k0.*alpha_G)./(alpha_L).*(alpha_G<1-1e-6);
    s1 = - s0.*alpha_G./(alpha_L) .*(alpha_G<1-1e-6);
    
    help1 = rho_L.*alpha_L.*k1+rho_G.*alpha_G.*k0;
    help2 = rho_L.*alpha_L.*s1+rho_G.*alpha_G.*s0;

    %  Averaging velocities.
    v_L = k1.*(I-help2)./help1+s1;
    v_G = k0.*(I-help2)./help1+s0;
    
%% Flux splitting
% indices
inL = 1:p.P-1;
inR = 2:p.P;
inF = 2:p.P;

% Compute approximate mixture sound velocity
% K_star =  (p.K-(p.K-1)*alpha_G);
K_star = 1;
den = alpha_G.*rho_L.*(1-K_star.*alpha_G);
omega = sqrt( P./max(1e-2,den) );
% omega = omega.*(omega<lambda_max) + lambda_max.*(omega>=lambda_max);

MachL = min(omega(inL),p.c_L)     .*(alpha_G(inL) <= .5) + ...
        min(omega(inL),p.c_G(inL)).*(alpha_G(inL) > .5);
MachR = min(omega(inR),p.c_L)     .*(alpha_G(inR) <= .5) + ...
        min(omega(inR),p.c_G(inR)).*(alpha_G(inR) > .5);
Mach = max(MachL,MachR);

% Velocity splitting liquid
VL_L = velSplit(v_L(inL), Mach, alpha_L(inR), 1);
VL_R = velSplit(v_L(inR), Mach, alpha_L(inL),-1);

% Velocity splitting gas
VG_L = velSplit(v_G(inL), Mach, alpha_G(inR), 1);
VG_R = velSplit(v_G(inR), Mach, alpha_G(inL),-1);

% Pressure splitting
vM_L = v_L(inL).*alpha_L(inL) + v_G(inL).*alpha_G(inL);
vM_R = v_L(inR).*alpha_L(inR) + v_G(inR).*alpha_G(inR);
PL = presSplit(vM_L,Mach,1);
PR = presSplit(vM_R,Mach,-1);

mL = alpha_L(inL).*rho_L(inL);
mR = alpha_L(inR).*rho_L(inR);
nL = alpha_G(inL).*rho_G(inL);
nR = alpha_G(inR).*rho_G(inR);
    
%% Domain Fluxes
% Liquid flux
flc(inF,1) = mL.*VL_L + mR.*VL_R;
flc(inF,2) = 0;
flc(inF,3) = mL.*VL_L.*v_L(inL) + mR.*VL_R.*v_L(inR);

% Gas Flux
fgc(inF,1)=0;
fgc(inF,2)= nL.*VG_L + nR.*VG_R;
fgc(inF,3)= nL.*VG_L.*v_G(inL) + nR.*VG_R.*v_G(inR);

% Pressure flux
fp(inF,1)= 0.0;
fp(inF,2)= 0.0;
fp(inF,3)= PL.*P(inL) + PR.*P(inR);    

F_G = p.g*(n+m).*sin(p.thetaVec) * N/Nf;
v_M = alpha_G.*v_G + alpha_L.*v_L;
switch p.frictionModel
    case 'turbulent'
        F_W = p.fVec.*(m+n).*abs(v_M).*v_M./(p.Dc-p.Dp);
    case 'Flatten'
        F_W = 32*v_M.*(alpha_G*p.eta_G+alpha_L*p.eta_L)./(p.Dc-p.Dp).^2;
    case 'linear'
        F_W = p.fVec.*(m+n).*v_M./(p.Dc-p.Dp);
    otherwise
        error('Friction Model not defined');
end
F_W = F_W*N/Nf;
%% Inlet fluxes
% Extrapolate to find the pressure at the inlet/bottom of the well.
p_bh = P(1)+0.5*(P(1)-P(2));
%     Compute massrates from the reservoir
[W_Gres,W_Lres] = LeftMassrates(p,p_bh);

flc(1,1)= (p.rho0_L*p.Q_Lbit + W_Lres)/p.A;
flc(1,2)= 0;
flc(1,3)= flc(1,1)*v_L(1);

fgc(1,1)= 0;
fgc(1,2)= (p.W_Gbit + W_Gres)/p.A;
fgc(1,3)= fgc(1,2)*v_G(1);

fp(1,1)= 0;
fp(1,2)= 0;     
% Extrapolate to find the pressure at the inlet/bottom of the well.
fp(1,3)= p_bh;

%% Outlet fluxes (open & closed conditions)
delta = 3e5;
% Choke Massrate
W_c = valveEquation(p,n(N),m(N),P(N),v_G(N),v_L(N),delta)/p.A; %Kg/s

if p.Z < eps % Closed choke
     flc(Nf,1) = 0;
     flc(Nf,2) = 0;
     flc(Nf,3) = 0;

     fgc(Nf,1) = 0;
     fgc(Nf,2) = 0;
     fgc(Nf,3) = 0;

     fp(Nf,1) = 0;
     fp(Nf,2) = 0;
     fp(Nf,3) = P(N) - 0.5*(P(N-1)-P(N));
elseif strcmp(p.valveModel,'constantP')   % || WHP
    flc(Nf,1) = m(N)*v_L(N);
    flc(Nf,2) = 0;
    flc(Nf,3) = flc(Nf,1)*v_L(N);

    fgc(Nf,1) = 0;
    fgc(Nf,2) = n(N)*v_G(N);
    fgc(Nf,3) = fgc(Nf,2)*v_G(N);

    fp(Nf,1) = 0;
    fp(Nf,2) = 0;
    fp(Nf,3) = p.p_s + 1*(P(N)-P(N-1));
else % Use choke model
    YG = .1*n(N).*v_G(N)/(m(N).*v_L(N)+n(N).*v_G(N));
    X_G =  YG+n(N).*v_G(N)/(m(N).*v_L(N)+n(N).*v_G(N));
    X_L = -YG+m(N).*v_L(N)/(m(N).*v_L(N)+n(N).*v_G(N));
    
    % Set fluxes
    flc(Nf,1) = X_L *W_c;
    flc(Nf,2) = 0;
    flc(Nf,3) = flc(Nf,1)*v_L(N);

    fgc(Nf,1) = 0;
    fgc(Nf,2) = X_G *W_c;
    fgc(Nf,3) = fgc(Nf,2)*v_G(N);

    fp(Nf,1) = 0;
    fp(Nf,2) = 0;

    fp(Nf,3) = P(N) - 0.5*(P(N-1)-P(N));

end  

%%
in = 1:p.P;
in1 = 2:p.P+1;

m = m - dtdx*( ...
         ( flc(in1,1) - flc(in,1) )...
        +( fgc(in1,1) - fgc(in,1) )...
        +( fp (in1,1) - fp (in,1) ));

    
n = n - dtdx*( ...
         ( flc(in1,2) - flc(in,2) )...
        +( fgc(in1,2) - fgc(in,2) )...
        +( fp (in1,2) - fp (in,2) ));

I = I - dtdx*( ...
         ( flc(in1,3) - flc(in,3) )...
        +( fgc(in1,3) - fgc(in,3) )...
        +( fp (in1,3) - fp (in,3) ))...
        -dt*( (F_W(in)) + F_G(in) );
     
     
end
uu = u2u(p,n,m,I);

end


function WHP = invValveEquation(p,n,m,v_G,v_L,rho_G,rho_L)

    if v_L <= 0
        WHP = p.p_s;
    else
        WHP = p.p_s + (p.A/p.Z/p.C_v)^2*( ...
            m*v_L/sqrt(rho_L) + n*v_G/(p.Y*sqrt(rho_G)) )^2;
    end

end




















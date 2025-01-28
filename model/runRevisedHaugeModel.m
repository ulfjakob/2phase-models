function [uu] = runRevisedHaugeModel(p, u0, T0)
% Parse states
L1 = u0(1); % [m] position of top of gas
Vg = u0(2); % [m3] volume of gas
Pc = u0(3); 
Pbh= u0(4);

% Local timestep
dt = .1;
Nt = ceil(p.dt/dt); % Number of local timesteps
dt = p.dt/Nt;

V = p.A*p.L;
MeanRho_G = 300;
MeanRho_L = (1e3+p.rho0_L)/2;
beta_L = p.c_L^2*p.rho0_L;
p.alpha_LS = 0.05;
F = p.f*1.1*10;
delta = 1e5; % [Pa]
C0 = 1/(1-p.alpha_LS);
alpha_d = .25;
gamma = 0.7;
t = T0;
Te = T0+p.dt;
k = 1;
while t < Te
    %% Resolve left BCs
    [W_Gres,W_Lres] = LeftMassrates(p,Pbh);
    rho_L = p.rho0_L + Pbh/p.c_L^2;       % Density of liquid
    rho_G = Pbh/p.c_G(1)^2;                % Density of gas
    W_L = W_Lres + p.Q_Lbit*p.rho0_L;
    W_G = W_Gres + p.W_Gbit;
    qL = W_L/rho_L;
    qG = W_G/rho_G;

    % Gas velocity
    vG = C0*(qL+qG)/p.A + p.v_inf;
    
    % Bubble length
    Lg = Vg/p.A/alpha_d;
    
    FG = MeanRho_L*F*(qG+qL)*abs(qG+qL)/p.A^2 + MeanRho_L*p.g;
    PL1 = Pc + L1*FG;
    
    % Convenience variables
    I1 = log( (PL1 + Lg*FG*(1-alpha_d))/( PL1 ) );
    barBeta = beta_L/( 1+p.A*beta_L/V*alpha_d/(gamma*FG)*I1 );
    Txe = p.A*vG/gamma*alpha_d*(1-C0*alpha_d)*I1;
    
    % Choke flow
    X_G = 1-1/p.Y*sqrt(Pc/p.c_G(end)^2/rho_L);
    C_K = 1-alpha_d*(L1<1)*(Vg>0)*C0*X_G;
    Qc =  p.C_v*p.Z/sqrt(MeanRho_L)*regRoot(Pc-p.p_s,delta) /C_K;
    
    %% Dynamics
    % Lumped pressure dynamics
    PcDot = barBeta/V*(qL + qG - Qc + Txe );
    VgDot = Txe - p.A*alpha_d*(1-C0*alpha_d)*I1/(gamma*FG*(1-alpha_d)) * PcDot ...
        + qG - alpha_d*(L1<1)*(Vg>0)*Qc;
    L1Dot = -C0*(qL+qG)/p.A - p.v_inf - 1/(p.A*alpha_d)*VgDot/2;
    
    Pc = Pc + dt*PcDot;
    Vg = max(0,Vg + dt*VgDot);
    L1 = max(0,L1 + dt*L1Dot);

    Pbh = Pc +FG*(p.L-Vg/p.A);
    
    % Update time
    t = t+dt;
    k = k+1;
end


%% Parse Output
uu(1) = L1; % [m] position of top of gas
uu(2) = Vg; % [m3] volume of gas
uu(3) = Pc;
uu(4) = Pbh;





























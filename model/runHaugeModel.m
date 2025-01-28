function [uu] = runHaugeModel(p, u0, T0)
% Parse states
h = u0(1); % [m] position of top of gas
V = u0(2); % [m3] volume of gas
Pc = u0(3); 
Pbh= u0(4);
Q_c = u0(5);

% Local timestep
dt = .1;
Nt = ceil(p.dt/dt); % Number of local timesteps
dt = p.dt/Nt;

tau=.2;
p.alphaL_dist = 0;
p.lambda = .2;
c_G = mean(p.c_G);
delta = 3e5; %[Pa]

t = T0;
Te = T0+p.dt;
k = 1;
while t < Te
    % Resolve left BCs
    [W_Gres,W_Lres] = LeftMassrates(p,Pbh);
    W_L = W_Lres + p.Q_Lbit*p.rho0_L;
    W_G = W_Gres + p.W_Gbit;
    Q_L = W_L/p.MeanRho_L;
    
    % Gas velocity
    vG = (Q_L)/p.A/(1-p.alpha_LS) + p.v_inf;
    
    % Choke flow
    Q_c_r = p.C_v*p.Z/sqrt(p.MeanRho_L) * regRoot( Pc-p.p_s, delta );
    Q_c = Q_c+(Q_c_r-Q_c)*tau;
    
    dV = Q_c-Q_L;
    vG_1 = vG + 1/(p.A*(1-p.alphaL_dist))* dV;

    h = max(0, h + dt*( - vG_1 ) );
    V = max(0, V + dt*dV );
    
    pg = p_idgas_param(V,p.mg,p);
    [Pc,Pbh] = find_pressures(pg,h,V,p.mg,Q_c,Q_L,p);
    
    % Update time
    t = t+dt;
    k = k+1;
end


%% Parse Output
uu(1) = h; % [m] position of top of gas
uu(2) = V; % [m3] volume of gas
uu(3) = Pc;
uu(4) = Pbh;
uu(5) = Q_c;





























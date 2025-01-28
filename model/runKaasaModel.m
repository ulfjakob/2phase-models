function uu = KaasaModel(p, u0, T0)

Pc  = u0(1);
Pbh  = u0(2);

% Local timestep
dt = .1;
Nt = ceil(p.dt/dt); % Number of local timesteps
dt = p.dt/Nt;

t = T0;
Te = T0+p.dt;

V = p.A*p.L;
c_G = mean(p.c_G);
MeanRho_L = (1e3+p.rho0_L)/2;
beta_L = p.c_L^2*p.rho0_L;
delta = 1e5; %[Pa]

k = 1;
while t < Te
    % Resolve left BCs
    [W_Gres,W_Lres] = LeftMassrates(p,Pbh);
    W_L = W_Lres + p.Q_Lbit*p.rho0_L;
    W_G = W_Gres + p.W_Gbit;
    Q_L = W_L/MeanRho_L;
    Q_G = W_G/Pbh*c_G^2;

    % Frictional pressure loss (Kaasa et al. 2008)
    dP_f = p.f*(W_L+W_G)^2 * 3e5;

    % Pressures
    Pbh = Pc +MeanRho_L*p.g*p.L + dP_f;

    % Choke flow
    Q_c = p.C_v*p.Z/sqrt(MeanRho_L) * regRoot( Pc-p.p_s, delta );

    Pc = Pc + beta_L/V*dt*( Q_L+Q_G - Q_c );


    
    % Update time
    t = t+dt;
    k = k+1;
end

%% Parse Output
uu(1) = Pc;
uu(2) = Pbh;






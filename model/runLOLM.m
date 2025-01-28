function uu = runLOLM(p, u0, T0)

m  = u0(1);
n  = u0(2);
% Pc = u0(3);
Pbh= u0(4);

% Local timestep
dt = .1;
Nt = ceil(p.dt/dt); % Number of local timesteps
dt = p.dt/Nt;

t = T0;
Te = T0+p.dt;

V = p.A*p.L;
c_G = mean(p.c_G);
delta = 1e5; %[Pa]
k = 1;
while t < Te
    % Compute choke pressure
% Trick to ensure numerical stability
    n = n.*(n>1e-8) + 1e-8.*(n<1e-8);
    m = m.*(m>1e-8);
      
    % Coefficients:
    a = 1/p.c_L.^2;
    b = p.rho0_L - m - c_G^2.*n*a;
    c = -p.rho0_L*c_G^2.*n;

    % Analytical solution:
    Pc = (-b+sqrt(b.^2 - 4*a.*c))/(2*a);  % Pressure

    %% Resolve left BCs
    [W_Gres,W_Lres] = LeftMassrates(p,Pbh);
    W_L = W_Lres + p.Q_Lbit*p.rho0_L;
    W_G = W_Gres + p.W_Gbit;

    % Frictional pressure loss (Kaasa et al. 2008)
    dP_f = p.f*(W_L+W_G)^2 * 3e5;

    % Pressures
    Pbh = Pc +(m+n)*p.g*p.L + dP_f + p.LOL_G;

    % Choke flow
    Wout = p.C_v*p.Z * sqrt(m+n*p.Y) * regRoot( Pc-p.p_s, delta );

    m = m + 1/V*dt*( W_L - m/(m+n) * Wout );
    n = n + 1/V*dt*( W_G - n/(m+n) * Wout );
    
    % Update time
    t = t+dt;
    k = k+1;
end


%% Parse Output
uu(1) = m;
uu(2) = n;
uu(3) = Pc;
uu(4) = Pbh;







function [uu,yy] = runRedDFM(p, uu0, T0)

% Number of cells
N = p.P;
uu = nan(N+2,1);

aG  = uu0(1:N);
Pc  = uu0(N+1);
Pbh = uu0(N+2);

% Local timestep
% ds = p.ds;
ds = p.L/N;
dt = min(ds/100,1);
Nt = ceil(p.dt/dt); % Number of local timesteps
dt_i = p.dt/Nt;

t = T0;
Te = T0+p.dt;

V = p.A*p.L;
beta_L = p.c_L^2*p.rho0_L;
C0 = 1/(1-p.alpha_LS);
F = 1;

delta = 1e5; % [Pa]
gamma = 1;

k = 1;
while t < Te
    %% Resolve left BCs
    [W_Gres,W_Lres] = LeftMassrates(p,Pbh);
    rho_L = p.rho0_L + Pbh/p.c_L^2;       % Density of liquid
    rho_G_bh = Pbh/p.c_G(1)^2;                % Density of gas
    W_L = W_Lres + p.Q_Lbit*p.rho0_L;
    W_G = W_Gres + p.W_Gbit;
    qL = W_L/rho_L;
    qG = W_G/rho_G_bh;
    
    %% Compute approximate pressure and velocity
    MeanRho_G = 30;
    rho_M = p.rho_L*(1-aG) + MeanRho_G*aG;
    
    % Frictional pressure loss
    QTOT = qG + qL;
    dFric_dx = F*p.fVec.* (QTOT/p.A).*abs((qG+qL)/p.A)./(p.Dc-p.Dp) .*rho_M;

    dGrav_dx = p.g*sin(p.thetaVec).*rho_M;
    dPdx = dFric_dx+dGrav_dx;
    P = Pc + flipud(cumtrapz(p.s,flipud(dPdx)));

    % Gas at boundary
    vG0 = 1/p.A/(1-p.alpha_LS)*(qL+qG) + p.v_inf;
    % Change in gas velocity by expansion
    I_vG = cumtrapz(p.s,aG*C0/gamma./max(P,1e4).* dPdx);
%     I_vG = cumsum(aG*C0/gamma./max(P,1e4).* dPdx *p.ds);
    vG  =  exp(I_vG).*(vG0);
    
    %% Lumped pressure dynamics
    if strcmp(p.valveModel,'Fjalstad')
        X_G = 1-1/p.Y*sqrt(Pc/p.c_G(end)^2/p.rho_L);
        C_K = max(1-aG(end)*X_G, 1e-3);
    %     C_K = 1;
        Qc =  p.C_v*p.Z/sqrt(p.rho_L)*regRoot(Pc-p.p_s,delta) /C_K;

        Txe = p.A*(vG(N)-vG0);
        barBeta = beta_L/(1+p.A*beta_L/V*sum(aG./max(gamma*P,1e4))*p.ds);
        % Lumped pressure dynamics
        PcDot = barBeta/V*( qL + qG - Qc + Txe );

        if Pc + dt_i*PcDot < p.p_s
            % Time step is too large: solve for and enforce steady state sol.
            Pc = p.p_s+regRootInv((qL+qG+Txe)*C_K/(p.C_v*p.Z)*sqrt(p.rho_L),delta);
        else
            Pc = Pc + dt_i*PcDot;
        end
        Pc = max(Pc,p.p_s);
    else 
        Pc = p.p_s;
        PcDot = 0;
    end
    
    Pbh = P(1);

    %% Source terms

    % Gas expansion source term
    E_G = aG.*(1-aG-p.alpha_LS)./((1-p.alpha_LS)*max(P,1e4)*gamma).*...
        (vG.*dPdx - PcDot);
    
    % PDE left BC.
    aG0 = qG/((qG+qL)/(1-p.alpha_LS)+p.A*p.v_inf);
    
    if abs(max(vG)*dt_i/p.ds)>1
        warning('CFL vioalation');
    end
    
    %% Dynamics
    in1 = 1:N;
%     % Second order Upwind
%     aG_pad = [2*aG0-aG(1); aG0; aG];
%     aGxp = 1/2*(3*aG_pad(in1+2)-4*aG_pad(in1+1)+aG_pad(in1));
%     % First order Upwind
    aG_pad = [aG0; aG];
    aGxp = aG_pad(in1+1)-aG_pad(in1);
    
    S = E_G;    % Void source term
    DaG = - dt_i/p.ds* vG.*aGxp  + dt_i* S;

    aG = aG + DaG;
    aG = max(aG,0);
    aG = min(aG,1);
    
    %%
    % Update time
    t = t+dt_i;
    k = k+1;
end

%% Parse Output
uu(1:N) = aG;
uu(N+1) = Pc;
uu(N+2) = Pbh;

yy.P = P;
yy.vG = vG;


















function p = parameters
%Generates parameters for the model
%%
p.Z      = 0.45;          % [-] Choke opening in [0,1]
p.C0     = 1.0;          % Drift-flux parameter
p.v_inf  = 0.00;
p.Y      = 0.13;          % Gas expansion factor used at choke
p.alpha_LS=1-1/p.C0;    % Drift-flux parameter
p.L      = 2530;         % Well Lenght
% p.Nz=p.L/p.dz;
p.R      = 8.314;            % [J/mol/K] Gas constant
p.Z_G    = 1;                % [-] Compressibility factor
p.T      = 273.15+12;              % [K] Temperature
%         p.SG     = 0.6;             % [-] Gas specific gravity
%         p.M      = p.SG*0.029;          % [kg/mol] Molar gas weight
p.M      = 0.0254;
p.R_G    = p.Z_G*p.R/p.M;    % [J/kg/K] Effective Specifig Gas constant
p.rho_L  = 1000;             % [kg/m^3] Liquid Density
p.g      = 9.81;             % [m/s^2] acceleration of gravity
p.theta  = pi/2;             % [rad] Well inclination
p.Dc     = 0.1524;              % [m] Casing diameter
p.Dp     = 0.0889;                % [m] Pipe outer diameter
p.A=pi*(p.Dc^2-p.Dp^2)/4;   % [m^2] cross sectional flow area
p.p_s    = 3e5;              % [Pa] downstram choke pressure
p.C_v   = 8.1e-3;           % [m^2] Choke constant
p.valveModel='Fjalstad';    % Valve Model: 'classical' or 'Fjalstad'
%         p.valveModel='constantP';    % Valve Model: 'classical' or 'Fjalstad'
p.slipLaw = 'Evje';
p.p_res  = 279e5;            % [Pa] Reservoir pressure
p.resType = 'PI';
p.k_L    = 0;            % [kg/s/Pa] Liquid production index
p.k_G    = 3.13e-7;           % [kg/s/Pa] Gas production index
p.P      = 50;
p.f      = 30.0e-3;           % Friction factor
p.fVec   = p.f*ones(1,p.P).';   %Space-dependent friction factor
p.frictionModel = 'turbulent'; %Friction model 'laminar' or 'turbulent'
p.W_Gbit = 0;  % [kg/s] Gas mas flow rate through bit
p.Q_Lbit =  13.33*1e-3;     % [m^3/s] Liquid volumetric flow rate through bit
p.DeltaT = 3600*2-10;
p.dt     = 10;       % [s] Simulation time step
p.t      = 1:p.dt:p.DeltaT; % [s] Simulation time vector
p.ds     = p.L/(p.P-1);  % [m] Spatial resolution of discretization
%         p.ds     = p.L/p.P;  % [m] Spatial resolution of discretization
p.s      = linspace(0,p.L,p.P); % [m] vector of spatial distributed nodes
p.thetaVec  = p.theta*ones(p.P,1); % Inclination
p.TVec=linspace(30,20,p.P).'+273.15; % [K] Temperature Profile
p.c_G=sqrt(p.R_G*p.TVec);    %Velocity of sound in gas
p.c_L=1000;                     %Velocitiy of sound in liquid
p.rho0_L=974.82;                  %[kg/m^3] Liquid density at vacum
p.K=1.0;                        %[-] Slip Law parameter (Evje law)
p.S=0.0;                        %[m/s] Slip Law parameter (Evje law)

%% Notes:
%
% Standard mass amounts, conversions:
% std m3 (0°C, 1 atm) = 4.461 58 E-02 kmol
% std m3 (15°C, 1 atm) = 4.229 32 E-02 kmol
% std ft3 (60°F, 1 atm) = 1.195 3 E-03 kmol
% Example:
% 1 std m3(15°C, 1 atm) of N2 = 0.0423 kmoles * 28 kg/kmoles = 1.18 kg
%
% Molar gas weight: when specifying molar gas weight, remember that we are
% after the molar gas weight of the molecule and not the atom. I.e. molar
% gas weight of N2 is 2*14.007 = 28.014
%
% Specific gravity gas: =approx 1.205kg/m3 =approx 29 g/mol @ standard
% conditions.
% Specific gravity of a gas is given by air under standard conditions.


















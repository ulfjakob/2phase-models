%This script runs the full model in open-loop. By default, the script generates
%a set of parameters by calling the script parameters.m. The model is
%usually initialized at the equilibrium, which is also computed in this
%script for the chosen set of parameters.

clear;
clc;

addpath model model\AUMSV' Model'\ model\SS' solvers'\
addpath Datasets\

% choose models to run
simModel = [1 2 3 4 5 6];
% simModel = [2];

%% Load state & parmaeters
p = parameters;
p.Y = 0.12;
p.C_v   = 9.1e-3;           % [m^2] Choke constant
p.LOL_G = 0;

load('kick2.mat');
O.s = linspace(0,p.L,166);

% Pressure drop contributions at SS 1-phas
dP_f = O.FF_fric(340,10)-O.FF_fric(175,10); % [bar]
dP_g = O.FF_grav(340,10)-O.FF_grav(175,10); % [bar]
dP_r = O.FF_rest(340,10)-O.FF_rest(175,10); % [bar]

p.p_res = rr.p_res(1)*1e5;
p.k_G = rr.k_G(1)*0;
% p.W_Gbit = 1;
p.Z = rr.Z(1);
p0 = O.BHCP(1);

% State Initialization
[qSS] = SSprofile(p,260e5);

%% Initialize state of other models
q10 = qSS;
[n,m,I] = u2u(p,q10);
%Outputs
[P, v_G, v_L, alpha_G, alpha_L, F_W, F_G, v_M,...
    rho_G, rho_L, Phi] =...
    variablesFromStates(p,n,m,I);
t = p.t;
Nt = length(t);
qq1 = [q10 zeros(3*p.P,Nt-1)];
q20 = [alpha_G; P(p.P); P(1)];
qq2 = [q20 zeros(p.P+2,Nt-1)];

ml = mean(m);
mg = mean(n);
q30 = [ml; mg; P(p.P); P(1)];
qq3 = [q30 zeros(4,Nt-1)];

q40 = [P(p.P);P(1)];
qq4 = [q40 zeros(2,Nt-1)];

h0 = p.L+100*1;
V0 = 1.0;
q50 = [h0; V0; P(p.P); P(1); v_M(end)*p.A];
qq5 = [q50 zeros(5,Nt-1)];
%
p.F = 750*p.L*6e1;
p.MeanRho_L = (1e3+p.rho0_L)/2;
p.mg = ( p.p_s + p.F*p.Q_Lbit + p.MeanRho_L*p.g*h0 )*...
    1e5*V0*p.M/(p.Z_G*p.R*p.T); % Mass of gas?
k_ks = 243;

h0 = p.L+100*0;
V0 = 0.0;
q60 = [h0; V0; P(p.P); P(1)];
qq6 = [q60 zeros(4,Nt-1)];


%% Plot Initialization

plotWhileRunning = 1;

% Figures
if 1
    hFig = figure(1);
    % Remove all child objects, but do not reset figure properties
    hAx(1) = subplot(211); title('WHP');   set(hAx(1), 'NextPlot', 'replacechildren');xlabel('Time[min]');
    hAx(2) = subplot(212); title('BHCP');    set(hAx(2), 'NextPlot', 'replacechildren'); xlabel('Time[min]')
    
    gFig = figure(2);
    % Remove all child objects, but do not reset figure properties
    hBx(1) = subplot(411); title('$Q_{GST}^L$');   set(hBx(1), 'NextPlot', 'replacechildren');xlabel('Time[min]');
    hBx(2) = subplot(412); title('HOLDUP');   set(hBx(2), 'NextPlot', 'replacechildren'); xlabel('Time[min]')
    ylim([0 1]);
    hBx(3) = subplot(413); title('$V_G^L,V_L^L$');  set(hBx(3), 'NextPlot', 'replacechildren'); xlabel('Time[min]')
    hBx(4) = subplot(414); title('$V_G^0,V_L^0$');  set(hBx(4), 'NextPlot', 'replacechildren'); xlabel('Time[min]')
end

%% Log init
Y1 = logInit(p);
Y2 = logInit(p);
Y3 = logInit(p);
Y4 = logInit(p);
Y5 = logInit(p);
Y6 = logInit(p);

%% Transient Simulation
tic
for k=1:Nt-1
    %% Interpolate controls and inputs
    p.p_res     = interp1(rr.t,rr.p_res*1e5,t(k));
    p.k_G       = interp1(rr.t,rr.k_G,t(k))*0;
    p.Z         = interp1(rr.t,rr.Z,t(k));
    
    p.p_s = O.p_sep(k);
    p.W_Gbit = O.GG_0(k)*1;
    
    %% Run sim and log variables for Sim1
    if  ismember(1,simModel)
        [qq1(:,k+1)] = runAUMSVsteps(p, qq1(:,k), t(k));

        %States
        [n,m,I] = u2u(p,qq1(:,k+1));
        %Outputs
        [P, v_G, v_L, alpha_G, alpha_L, F_W, F_G, v_M,...
            rho_G, rho_L, Phi] =...
            variablesFromStates(p,n,m,I);

        W_G_c = n(p.P)*v_G(p.P)*p.A;
        Y1.GG_L(k) = W_G_c;
        Y1.WHP(k) = P(p.P);
        Y1.BHCP(k) = P(1)+0*(P(2)-P(1));
        Y1.VG_L(k) = v_G(p.P);
        Y1.VL_L(k) = v_L(p.P);
        Y1.VG_0(k) = v_G(1);
        Y1.VL_0(k) = v_L(1);

        Y1.HOLDUP(:,k) = 1-alpha_G;
    end
    
    %% Run sim and log variables for Sim2
    if  ismember(2,simModel)
        [qq2(:,k+1)] = runRedDFM(p, qq2(:,k), t(k));
        Y2.WHP(k)       = qq2(p.P+1,k);
        Y2.BHCP(k)      = qq2(p.P+2,k);
        Y2.HOLDUP(:,k)  = 1-qq2(1:p.P,k);
    end
    
    %% Run sim and log variables for Sim3
    if  ismember(3,simModel)
        [qq3(:,k+1)] = runLOLM(p, qq3(:,k), t(k));
        Y3.WHP(k)    = qq3(3,k);
        Y3.BHCP(k)   = qq3(4,k);
    end
    
    %% Run sim and log variables for Sim4
    if  ismember(4,simModel)
        [qq4(:,k+1)] = runKaasaModel(p, qq4(:,k), t(k));
        Y4.WHP(k)       = qq4(1,k);
        Y4.BHCP(k)      = qq4(2,k);
    end
    
    %% Run sim and log variables for Sim5
    if k == k_ks
        qq5(:,k) = q50;
    end
    if  ismember(5,simModel)  && k >= k_ks
        [qq5(:,k+1)] = runHaugeModel(p, qq5(:,k), t(k));
        Y5.WHP(k)       = qq5(3,k);
        Y5.BHCP(k)      = qq5(4,k);
%         qq5(1:2,k)
    end
    
    %% Run sim and log variables for Sim6
    if k < k_ks
        qq6(1,k) = p.L;
    end
    if  ismember(6,simModel)
        [qq6(:,k+1)] = runRevisedHaugeModel(p, qq6(:,k), t(k));
        Y6.WHP(k)       = qq6(3,k);
        Y6.BHCP(k)      = qq6(4,k);
    end
    
    %% 
    if plotWhileRunning && mod(k,10) == 0
        Tnom = 60;
        plot(hAx(1), t/Tnom, O.WHP(1:Nt)*1e-5, t/Tnom, Y1.WHP*1e-5, ...
            t/Tnom, Y2.WHP*1e-5, t/Tnom, Y3.WHP*1e-5, ...
            t/Tnom, Y4.WHP*1e-5, t/Tnom, Y5.WHP*1e-5, ...
            t/Tnom, Y6.WHP*1e-5);
        plot(hAx(2), t/Tnom, O.BHCP(1:Nt)*1e-5, t/Tnom, Y1.BHCP*1e-5, ...
            t/Tnom, Y2.BHCP*1e-5, t/Tnom, Y3.BHCP*1e-5, ...
            t/Tnom, Y4.BHCP*1e-5,t/Tnom, Y5.BHCP*1e-5, ...
            t/Tnom, Y6.BHCP*1e-5);

        plot(hBx(1), t/Tnom, O.GG_L(1:Nt), t/Tnom, Y1.GG_L);
        plot(hBx(2), O.s, O.HOLDUP(175:340,k), p.s, Y1.HOLDUP(:,k), ...
            p.s, Y2.HOLDUP(:,k));
        plot(hBx(3), t/Tnom, O.VG_L(1:Nt), t/Tnom, Y1.VG_L,'--',...
            t/Tnom, O.VL_L(1:Nt), t/Tnom, Y1.VL_L,'--');
        plot(hBx(4), t/Tnom, O.VG_0(1:Nt), t/Tnom, Y1.VG_0,'--',...
            t/Tnom, O.VL_0(1:Nt), t/Tnom, Y1.VL_0,'--');
        
        drawnow;
    end
    
    
end
toc

%% Save
% uFinal=q(:,k);
% save output/uFinal uFinal
% save output/lastRun u t

%%
plot_marker = nan(size(t));
plot_marker(k_ks:k_ks+1) = [0 1e3];
plot_marker(50*6:50*6+1) = [0 1e3];

clc
lineStyles = {'-', '--', '-.', '-.'};
% myColors = lines(4);
myColors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

figure(5); clf;
subplot(211)
t1 = 35;
t2 = 100;

plot(t/Tnom, O.WHP(1:Nt)*1e-5, 'linestyle','-' ,'color',myColors(1,:),'linewidth',2);
hold on;
% plot(t/Tnom, plot_marker,':k');
legend_str = {'OLGA'};
ylim([0 80]);
xlim([t1 t2])
ylabel('WHP [bar]');
xlabel('Time [min]');
legend(legend_str,'location','NorthWest');
% columnlegend(2,legend_str,'NorthWest');

subplot(212)
plot(t/Tnom, O.BHCP(1:Nt)*1e-5, 'linestyle','-' ,'color',myColors(1,:),'linewidth',2);
hold on;
% plot(t/Tnom, plot_marker,':k');
ylim([230 320]);
xlim([t1 t2])
ylabel('BHCP [bar]');
xlabel('Time [min]');

%
figure(5); 
% clf;
subplot(211)
t1 = 35;
t2 = 100;

plot(t/Tnom, Y1.WHP(1:Nt)*1e-5,'linestyle','--','color',myColors(2,:),'linewidth',2);
hold on;
plot(t/Tnom, Y2.WHP(1:Nt)*1e-5,'linestyle','-.' ,'color',myColors(3,:),'linewidth',2);
plot(t/Tnom, plot_marker,':k');
legend_str = {'OLGA','Mec DFM','Red DFM'};
ylim([0 80]);
xlim([t1 t2])
ylabel('WHP [bar]');
xlabel('Time [min]');
legend(legend_str,'location','NorthWest');

subplot(212)
plot(t/Tnom, Y1.BHCP(1:Nt)*1e-5,'linestyle','--','color',myColors(2,:),'linewidth',2);
hold on;
plot(t/Tnom, Y2.BHCP(1:Nt)*1e-5,'linestyle','-.' ,'color',myColors(3,:),'linewidth',2);
plot(t/Tnom, plot_marker,':k');
ylim([230 320]);
xlim([t1 t2])
ylabel('BHCP [bar]');
xlabel('Time [min]');


%
figure(7); clf;
subplot(211)
t1 = 35;
t2 = 100;

plot(t/Tnom, Y3.WHP(1:Nt)*1e-5,'linestyle','-','color',myColors(4,:),'linewidth',2);
hold on;
plot(t/Tnom, Y4.WHP(1:Nt)*1e-5,'linestyle','--','color',myColors(5,:),'linewidth',2);
plot(t/Tnom, Y5.WHP(1:Nt)*1e-5,'linestyle','-.','color',myColors(6,:),'linewidth',2);
% plot(t/Tnom, Y6.WHP(1:Nt)*1e-5,'linestyle',':','color',myColors(7,:),'linewidth',2);
plot(t/Tnom, plot_marker,':k');
% legend_str = {'LOL mod','1-ph mod','Hauge','Rev. Hauge'};
legend_str = {'LOL mod','1-ph mod','Hauge'};
ylim([0 80]);
xlim([t1 t2])
ylabel('WHP [bar]');
xlabel('Time [min]');
legend(legend_str,'location','NorthWest');

subplot(212)
plot(t/Tnom, Y3.BHCP(1:Nt)*1e-5,'linestyle','-','color',myColors(4,:),'linewidth',2);
hold on;
plot(t/Tnom, Y4.BHCP(1:Nt)*1e-5,'linestyle','--','color',myColors(5,:),'linewidth',2);
plot(t/Tnom, Y5.BHCP(1:Nt)*1e-5,'linestyle','-.','color',myColors(6,:),'linewidth',2);
% plot(t/Tnom, Y6.BHCP(1:Nt)*1e-5,'linestyle',':','color',myColors(7,:),'linewidth',2);
plot(t/Tnom, plot_marker,':k');
ylim([200 320]);
xlim([t1 t2])
ylabel('BHCP [bar]');
xlabel('Time [min]');


%%
% gcf = figure(5);
% pdfmatlabfrag2(gcf, 'PressureTrend_kick_HiFiDFM.pdf');
% gcf = figure(6);
% pdfmatlabfrag2(gcf, 'PressureTrend_kick_DFM.pdf');
% gcf = figure(7);
% pdfmatlabfrag2(gcf, 'PressureTrend_kick_LOL.pdf');














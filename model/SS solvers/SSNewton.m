function [uu,RelError] = SSNewton(p,uu,itMax)
% function uu = SSNewton(p,uu)
%
% Quick, efficient way to find a steady-state.
% 
% Syntax 
% 
% [uu] = SSNewton(p,uu) finds a steady-state of the parameter struct p by
% using newton iterations with an initial guess uu, consisting of the
% conservative states vertically concatenated.
%
% [uu] = SSNewton(p,uu,itMax) does so with a maximum of itMax iterations
% 
% [uu,RelError] = SSNewton(p,uu,itMax) also returns a weighted norm of
% the left hand side of the steady-state equations, which is a measure of
% how close to the steady-state the algorithm has converged. 
% 
% Description
% 
% This function simply consists of running the transient model, with an
% infinite time-step, which corresponds to the discretized steady-state 
% equations. It has the advantage of computing a steady-state solution
% which is consistent with the transient model. However, it is rigorously
% less accurate than solving the steady-state equations with a Matlab ode45
% solver. 

%   Ulf Jakob F. Aarsnes    01.04.2014

if nargin == 2
    itMax = 5;
end

% Find SS Jacobian by setting timestep to inf hence removing the effect
% of the time derivative term.
dt = p.dt;
p.dt = inf;
if isfield(p,'modelType') && strcmp(p.modelType,'Riemann')
    [f,J,RelError]=fDFModel_without_reduction(0,uu,uu,p);
else
    [f,J,RelError]=fDFModel_full(0,uu,uu,p);
end
p.dt = dt;

% Use numerical jacobian if analytics is undefined
if isnan(rcond(J)) || rcond(J) < 1e-16
    if isfield(p,'modelType') && strcmp(p.modelType,'Riemann')
        [f,~,RelError]=fDFModel_without_reduction(0,uu,uu,p);
        J = numjac(@(tt,x)fDFModel_without_reduction(tt,x,x,p),0,uu,f,1,[],'off');
    else
        [f,~,RelError]=fDFModel_full(0,uu,uu,p);
        J = numjac(@(tt,x)fDFModel_full(tt,x,x,p),0,uu,f,1,[],'off');
    end
end

it = 0;
while RelError > 1e-4 && it < itMax
   
    uu = uu - J\f;
    %% Project the states into physically meaningfull values. 
    % Doesn't seem to have significant impact on results but sometimes
    % keeps the simulator from crashing.
    uu(1:2*p.P) = max(uu(1:2*p.P),0);
    uu = real(uu);
        
    it = it+1;
    
    % Find SS Jacobian by setting timestep to inf hence removing the effect
    % of the time derivative term.
    dt = p.dt;
    p.dt = inf;
    if isfield(p,'modelType') && strcmp(p.modelType,'Riemann')
        [f,J,RelError]=fDFModel_without_reduction(0,uu,uu,p);
    else
        [f,J,RelError]=fDFModel_full(0,uu,uu,p);
    end
    p.dt = dt;
    if isfield(p,'modelType') && strcmp(p.modelType,'Riemann')
        if isnan(rcond(J)) || rcond(J) < 1e-16
            [f,~,RelError]=fDFModel_without_reduction(0,uu,uu,p);
            J = numjac(@(tt,x)fDFModel_without_reduction(tt,x,x,p),0,uu,f,1,[],'off');
        end
    else
        if isnan(rcond(J)) || rcond(J) < 1e-16
            [f,~,RelError]=fDFModel_full(0,uu,uu,p);
            J = numjac(@(tt,x)fDFModel_full(tt,x,x,p),0,uu,f,1,[],'off');
        end
    end


end



















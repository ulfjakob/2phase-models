function [P0SS,qSS, flag] = SSfinder(p,n_SS,dOc)
%% [P0SS,uSS flag] = SSfinder(p,n_SS,dOc)
% Returns steady states (SS) of the parameter struct p.
% Inputs: 
%       p:      well parameters struct.
%       n_SS:   [Optional] which SS' are should be found.
%           n_SS = [1] : Return all SS'
%           n_SS = [2] : Return stable UB SS if exists. Else OB.
%       dOc:    [optional] Discretized or dist solution?
%           dOc: = [1] : distributed
%           dOc: = [2] : discretized
%
% Outputs: 
%       p0SS:       [Pa] Vector of p0 corresponding to the SS'.
%           p0SS(1):    UB stable (NaN if non-existent).
%           p0SS(2):    UB unstable (NaN if non-existent).
%           p0SS(3):    overbalanced (NaN if non-existent).
%       uSS:        SS States vector.
%           uSS(:,1):    UB stable (NaN if non-existent).
%           uSS(:,2):    UB unstable (NaN if non-existent).
%           uSS(:,3):    overbalanced (NaN if non-existent).
%       flag:       [text] Info about SS, if any.
%
% Ulf Jakob F. Aarnes, 14-05-2013
% Updated Flo Di Meglio 31/01/2013
% Updated Ulf Jakob F. Aarsnes 03/04/2014
%% Parse input
if nargin == 0;
    p = parameters;
    n_SS = 1;
    dOc = 1;
elseif nargin == 1;
    n_SS = 1;
    dOc = 1;
elseif nargin == 2
    dOc = 1;
end

%% Init output
flag = [];
P0SS = nan*zeros(1,3);
qSS = nan*zeros(3*p.P,3);

%% Check h3 at reservoir pressure
h3 = H3(p,p.p_res);

% h3<0 means 1 UB SS
if h3 < 0 || n_SS == 2
    for p0 = [.8 .55 .95]*p.p_res
        try
            [P0SS(1),qSS(:,1)] = getP0uSS(p,p0);
            break;
        catch err
            if (strcmp(err.identifier,'cmpSS:PTendsTo0') || ...
                    strcmp(err.identifier,'full:imagNan') )
                % All is good
            else
                rethrow(err);
            end
        end
    end
end
if h3 < 0 || ~isnan(P0SS(1))
    return;
end
%% 1 overbalanced SS (+ possibly 2 UB).
% Find overbalanced SS
[P0SS(3),qSS(:,3)] = getP0uSS(p,p.p_res+1e5);

%% Try to find negative H3
[p0,h3] = fminbnd(@(x)H3(p,x),p.p_res/2,p.p_res);

if h3 > 0
    % No UB SS'.
    return;
end

% Find  Unstable UB SS
p0_lb = p0;
p0_ub = p.p_res;
options = optimset('TolX',1e-6);
[P0SS(2),fval] = fzero(@(x)H3(p,x),[p0_lb,p0_ub],options);
if dOc == 1; [qSS(:,2)] = SSprofile(p,P0SS(2));
else [P0SS(2),qSS(:,2)] = getP0uSS(p,P0SS(2));
end

% Find  stable UB SS
for p0_lb = [.8 .7 .9 .6 .95 .5]*p0
    if inf > H3(p,p0_lb) > 0
        break
    end
end
p0_ub = p0;
[P0SS(1),fval] = fzero(@(x)H3(p,x),[p0_lb,p0_ub],options);
if dOc == 1; [qSS(:,1)] = SSprofile(p,P0SS(1));
else [P0SS(1),qSS(:,1)] = getP0uSS(p,P0SS(1));
end

end

function [p0,uSS] = getP0uSS(p,p0)
    [uSS] = SSprofile(p,p0);
    [uSS,RelError] = SSNewton(p,uSS,10);
    if RelError > 1e16 || isnan(RelError)
        p0 = nan;
        uSS = nan*zeros(3*p.P,1);
        return
    end
    [n,m,I] = u2u(p,uSS);
    P = variablesFromStates(p,n,m,I);
    p0 = P(1);
end

function h3 = H3(p,p0)
    try
        [uSS] = SSprofile(p,p0);
    catch err
        if (strcmp(err.identifier,'cmpSS:PTendsTo0'))
            h3 = inf;
            return
        else
            rethrow(err);
        end
    end
    dt = p.dt;
    p.dt = inf;
    f=fDFModel_full(0,uSS,uSS,p);
    p.dt = dt;
    h3 = f(end); 
end

























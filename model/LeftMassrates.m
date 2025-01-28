function [W_Gres,W_Lres,dW_GresdP,dW_LresdP] = ...
    LeftMassrates(p,P1)

%% LeftMassrates
% 
% Computes the reservoir inflow mass rates of gas, liquid, and their 
% derivatives
% 
% Syntax
% 
%  [W_Gres,W_Lres,dW_GresdP,dW_LresdP] = ...
%     LeftMassrates(p,P1) computes the following quantities:
% 
%  - W_Gres : gas mass flow rate from the reservoir [kg/s]
%  - W_Gres : liquid mass flow rate from the reservoir [kg/s]
%  - dW_GresdP : partial derivative of W_Gres w.r.t. Bottom hole pressure P
%  - dW_GresdP : partial derivative of W_Lres w.r.t. Bottom hole pressure P
% 
% for a given parameter struct p and a Bottom Hole Pressure P1. 
% 
% WARNING: the computation depends on the value of the parameter p.restype,
% which can take the value 'table', in the case of an IPR reservoir model
% or anything else for a Productivity Index type of model. In the latter
% case, the inflow coefficients p.k_G and p.k_L need to be defined in the
% parameter struct. 
% 
% 

if strcmpi(p.resType,'table')
    [W_Gres dW_GresdP] = TABLE_IPR(p.p_res-P1, p.tableIPR);
    W_Lres = 0;
    dW_LresdP = 0;
else
    % Regularized Pressure Drawdown and Jacobian
    [DP_res, dDP_resdP] = regMax(p.p_res-P1,.5e5);
    W_Gres = p.k_G.*DP_res;
    dW_GresdP = p.k_G*dDP_resdP;
    W_Lres = p.k_L.*DP_res;
    dW_LresdP = p.k_L*dDP_resdP;
end


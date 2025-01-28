function [y,dydx] = regRoot(x,delta) 
%Anti-symmetric square root approximation with finite derivative in the
%origin
%
% Inputs:
%       x
%       delta:  Range of significant deviation from sqrt(abs(x))*sgn(x), 
%               default is 0.01. 
%   output:
%       y
%       dy
%   algorithm:
%       y := x/(x*x+delta*delta)^0.25;
%       dydx := 0.5*(x.*x+2*delta*delta)./((x.*x+delta*delta).^1.25);
% url: 
%       http://www.maplesoft.com/documentation_center/online_manuals/modelica/Modelica_Fluid_Utilities.html#Modelica.Fluid.Utilities.regRoot

y = x./(x.*x+delta*delta).^0.25;
dydx = 0.5*(x.*x+2*delta*delta)./((x.*x+delta*delta).^1.25);

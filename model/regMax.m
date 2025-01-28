%% function [f,J]=regMax(x,k)
% Returns regularized approximation of max and it's derivative.
% The approximatino is valid when x >> k.
% f =~ max(x,0)
%
% Mar-2014-U.J.F. Aarsnes: Created based on Modelicas regRoot


function [f,J]=regMax(x,k)

f = max(0,x.^3./(x.^2+k^2));
x = (x>0).*x;
J = max((3*k^2*x.^2+x.^4)./(k^2+x.^2).^2,0);



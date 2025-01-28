function [x] = regRootInv(y,delta) 
% Inverse of the regRoot function
%
x = sign(y).*sqrt( 1/2*abs((y.^4 + sqrt(abs(y.^8+4*delta^2*y.^4)) )) );

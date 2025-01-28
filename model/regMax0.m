function [f,J]=regMax0(x,k)
f = max(x,k);
J = (x>k);
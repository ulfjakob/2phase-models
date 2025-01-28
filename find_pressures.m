function [ pc,pbit ] = find_pressures( pg,L1,V,mg,qc,qp,p )
%FIND_PRESSURES Calculate the pressure at choke and bit.
%   The choke and bit pressure, pc and pbit [bar]

hg=2*V/(p.A*(1-p.alphaL_dist)); % [m] initial height of gas

pfricG=p.F*qc*hg/p.L;
pfric1=p.F*qc*L1/p.L;
L2=p.L-hg-L1;
pfric2=p.F*qp*L2/p.L;

pc = pg - p.MeanRho_L*p.g*L1 - p.lambda*pfricG - pfric1 - p.lambda*p.g*(1+p.alphaL_dist)/...
    (1-p.alphaL_dist)*p.MeanRho_L*V/p.A + mg*1e-5/p.A;
pbit = pg + p.MeanRho_L*p.g*L2+(1-p.lambda)*pfricG - pfric2 + (1-p.lambda)*p.g*...
    (1+p.alphaL_dist)/(1-p.alphaL_dist)*p.MeanRho_L*V/p.A + mg*1e-5/p.A;

end


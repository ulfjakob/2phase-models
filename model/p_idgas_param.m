function [ pgas ] = p_idgas_param(Vgas,mg,p)
%P_IDGAS_PARAM computes the pressure of a gas bubble using the ideal gas law.
%   The pressure of the gas, pgas [bar] is computed from the factor Z
%   [0-1], the mass of gas, mg [kg], the molar mass of the gas, Mg
%   [kg/mol], the temperature, T [K], and the volume of gas, Vgas [m3]. p
%   is a parameter struct.

pgas=p.Z_G*mg/p.M*p.R*p.T./Vgas*1e-5; % [bar] pressure with ideal gas law

end


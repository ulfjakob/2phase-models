function Vs = velSplit(v,Mach,alpha,sign)
    if sign ~= 1 && sign ~= -1
        error('sign must be 1 or -1');
    end
    oVs = sign*alpha.*(v+sign*Mach).^2./(4*Mach) + (1-alpha).*(v+sign*abs(v))/2;
    lVs = .5*(v+sign*abs(v));
    
    % Choose velocity depending on whether we exceed Mach or not.
    Vs =  oVs .* (abs(v) <= Mach) ...
        + lVs .* (abs(v) >  Mach);
end




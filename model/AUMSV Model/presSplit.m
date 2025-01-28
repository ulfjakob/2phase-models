function Ps = presSplit(v,Mach,sign)
    if sign ~= 1 && sign ~= -1
        error('sign must be 1 or -1');
    end
    oPs = sign*(v+sign*Mach).^2./(4*Mach).*(sign*2-v./Mach)./Mach;
    lPs = (v+sign*abs(v))./(2*v);

    % Choose velocity depending on whether we exceed Mach or not.
    Ps =  oPs .* (abs(v) <= Mach) ...
        + lPs .* (abs(v) >  Mach);
end

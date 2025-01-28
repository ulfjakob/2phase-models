function Y = logInit(p)

Y.WHP = nan(size(p.t));
Y.BHCP = nan(size(p.t));
Y.GG_L = nan(size(p.t));
Y.VG_L = nan(size(p.t));
Y.VL_L = nan(size(p.t));
Y.VG_0 = nan(size(p.t));
Y.VL_0 = nan(size(p.t));

% Profile
Y.HOLDUP = nan(p.P,numel(p.t));
Y.P      = nan(p.P,numel(p.t));
Y.vG     = nan(p.P,numel(p.t));
end
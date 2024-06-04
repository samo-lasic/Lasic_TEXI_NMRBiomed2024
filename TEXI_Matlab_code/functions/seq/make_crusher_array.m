function [qc_array, gc_array, dz] = make_crusher_array(Nc, dfds, delta_c, dz_min, dz_max)
% dz_max only used when Nc > 2
if nargin < 5
    dz_max = 10 * 1e-3;
end

qc_max = pi*1e4 / (dz_min * 1e3); % scale_dz relative to 1 mm slice

%dfds = 6; % ﻿Lasič et al., Magn Reson Med. 2018;79(4):2228–35.
dz_min = butterfly2slice(qc_max, dfds) * 1;

%dz = [dz_min*linspace(1,40,6) Inf];
%dz = [dz_min*logspace(0,log10(10),10) Inf];

if Nc > 2
    dz = dz_min+(dz_max-dz_min)*linspace(0,1,Nc-1).^2;
else
    dz = dz_min;
end
dz = [dz Inf];

[qc_array, gc_array] = slice2butterfly(dz, delta_c, dfds);
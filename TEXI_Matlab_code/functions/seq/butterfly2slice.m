function dz = butterfly2slice(q_butterfly, BW_delta_slice)
% q_min_butterfly = 4*pi/dz+q_slice - a combination of crusher and slice gradients
% q_slice = gmr*g_slice*delta_slice/2 - dephasing from slice gradient

% 2*pi*rf_BW = gmr*g_slice*dz - RF bandwidth and slice thickness
%BW_delta_slice = 6; % ﻿Lasič S, Lundell H, Topgaard D, Dyrby TB. Effects of imaging gradients in sequences with varying longitudinal storage time—Case of diffusion exchange imaging. Magn Reson Med. 2018;79(4):2228–35.

dz = pi*(4+BW_delta_slice)./q_butterfly;
end
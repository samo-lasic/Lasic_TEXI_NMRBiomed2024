
function [q_min_butterfly, g_min_butterfly] = slice2butterfly(dz, delta_butterfly, BW_180_delta_slice)
% q_min_butterfly = 4*pi/dz+q_slice - a combination of crusher and slice gradients
% q_slice = gmr*g_slice*delta_slice/2 - dephasing from slice gradient
% 2*pi*BW_180 = gmr*g_slice*dz - RF bandwidth and slice thickness
%BW_180_delta_slice = 6; % ﻿Lasič S, Lundell H, Topgaard D, Dyrby TB. Effects of imaging gradients in sequences with varying longitudinal storage time—Case of diffusion exchange imaging. Magn Reson Med. 2018;79(4):2228–35.
gmr = 26.75e7;
q_min_butterfly = pi*(4+BW_180_delta_slice)./dz;
g_min_butterfly = q_min_butterfly/gmr/delta_butterfly;
end

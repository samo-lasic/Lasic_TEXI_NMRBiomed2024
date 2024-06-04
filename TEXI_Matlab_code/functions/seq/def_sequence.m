function [g, q] = def_sequence(delta, delta_max, G1, G2, tm, gc, delta_c, delta_pi, dt)
% allows for different transverse encoding blocks 

% if nargin < 10
%     polarity = 1;
% elseif invert_polarity == 1
%     polarity = -1;
% else
%     polarity = 1;
% end

gmr = 26.75e7;
gc_pulse = gc*ones(round(delta_c/dt),1);
g180 = zeros(round(delta_pi/dt),1);

gT = ones(round(delta/dt),1);
g0 = zeros(round((delta_max-delta)/dt),1);
gT = [gT; g180; -gT; g0];

gmix = [gc_pulse; zeros(round(tm/dt),1); -gc_pulse];

g = [G1 * gT; gmix; -G2 * flipud(gT)];

q = gmr*cumsum(g)*dt;

end
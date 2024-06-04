function [g, q] = def_FEXI(delta, delta_max, G1, G2, tm, gc, delta_c, delta_pi, dt)
gmr = 26.75e7;
gc_pulse = gc*ones(round(delta_c/dt),1);
g180 = zeros(round(delta_pi/dt),1);

gT = ones(round(delta/dt),1);
g0 = zeros(round((delta_max-delta)/dt),1);
gT1 = G1*[gT; g180; -gT; g0];
gT2 = G2*[gT; g180; -gT; g0];

gmix = [gc_pulse; zeros(round(tm/dt),1); -gc_pulse];

g = [gT1; gmix; -flipud(gT2)];
%gn = [gT; gmix; flipud(gT)];

q = gmr*cumsum(g)*dt;
%qn = gmr*cumsum(gn)*dt;

end
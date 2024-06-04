function s = Ning_h3_norm_fit2data(m, xps)
%Ning Ning, Nilsson M, Lasič S, Westin C-F, Rathi Y. Cumulant expansions for measuring water exchange using diffusion MRI. J Chem Phys. 2018;148(7).
% modification of Ning_h by adding one more cumulant

MD  = m(1);
V   = m(2);
k   = m(3);
c3  = m(4);
% S0  = m(5:end);

b = xps.b;
b2 = xps.b.^2;
b3 = xps.b.^3;
h = 2 * sum(exp(-k*xps.t).*xps.q4) * xps.dt ./ b2;
s = exp(-b.*MD + (b2/2 * V - b3/6 * c3).*h);

end
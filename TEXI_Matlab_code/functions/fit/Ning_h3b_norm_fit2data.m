function s = Ning_h3b_norm_fit2data(m, xps)
%Ning Ning, Nilsson M, Lasič S, Westin C-F, Rathi Y. Cumulant expansions for measuring water exchange using diffusion MRI. J Chem Phys. 2018;148(7).
% modification of Ning_h3 by scaling exchange sensitivity h with b-value

MD  = m(1);
V   = m(2);
k   = m(3);
c3  = m(4);
% S0  = m(5:end);

h2 = sum(exp(-k*xps.t).*xps.q4) * xps.dt;
h3 = h2.^(3/2);
s = exp(-xps.b.*MD + (h2/2 * V - h3/6 * c3));

end
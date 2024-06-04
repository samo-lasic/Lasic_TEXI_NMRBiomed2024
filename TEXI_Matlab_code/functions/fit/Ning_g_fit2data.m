function s = Ning_g_fit2data(m, xps)
%Ning Ning, Nilsson M, Lasič S, Westin C-F, Rathi Y. Cumulant expansions for measuring water exchange using diffusion MRI. J Chem Phys. 2018;148(7).
% modification of Ning_h by using the gamma distribution function

MD  = m(1);
V   = m(2);
k   = m(3);
S0  = m(4:end);

h = 2 * sum(exp(-k*xps.t).*xps.q4) * xps.dt ./ xps.b.^2;

Vh = V * h;
s = (1 + xps.b .* Vh/MD).^(-MD^2 ./ Vh);

s = s .* S0(xps.tm_ind);
end
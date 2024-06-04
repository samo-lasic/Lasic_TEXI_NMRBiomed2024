function s = Ning_Gamma_fit2data(m, xps)
%Ning Ning, Nilsson M, Lasič S, Westin C-F, Rathi Y. Cumulant expansions for measuring water exchange using diffusion MRI. J Chem Phys. 2018;148(7).
% using exchange weighting time

MD  = m(1);
V   = m(2);
k   = m(3);
S0  = m(4:end);

s = exp( -xps.b.*MD + xps.b.^2/2 .* (1 - k*xps.Gamma)*V );

s = s .* S0(xps.tm_ind);
end
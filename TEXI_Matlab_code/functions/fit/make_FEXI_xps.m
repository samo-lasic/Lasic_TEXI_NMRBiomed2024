function xps = make_FEXI_xps(seq);
% xps.n            = 6;
% xps.mde_b1       = [0 0  1 1  1 1]' * 0.9e9;
% xps.mde_b2       = [0 1  0 1  0 1]' * 1.0e9;
% xps.mde_tm12     = [0 0  0 0  1 1]' * 0.3;
% xps.s_ind        = [1 1  2 2  3 3]';
% xps.mde_tm12_ind = [1 1  1 1  2 2]';
% xps.mde_b1_ind   = [1 1  2 2  2 2]';
% xps.mde_b2_ind   = [1 2  1 2  1 2]';

xps.n = seq.n;
xps.mde_b1 = seq.b1';
xps.mde_b2 = seq.b2';
xps.mde_tm12 = seq.tm';

xps.mde_tm12_ind = ones(1,seq.n);
u = unique(seq.tm);
for c = 1:length(u)
    ind = find(seq.tm == u(c));
    xps.mde_tm12_ind(1,ind) = c;
end

xps.s_ind = seq.s_ind;

xps.mde_b1_ind   = ones(1,xps.n);
xps.mde_b2_ind   = ones(1,xps.n);

u = unique(seq.b1);
for c = 1:length(u)
    ind = find(seq.b1 == u(c));
    xps.mde_b1_ind(ind) = c;
end

u = unique(seq.b2);
for c = 1:length(u)
    ind = find(seq.b2 == u(c));
    xps.mde_b2_ind(ind) = c;
end

xps.tm_ind = xps.mde_tm12_ind;

end

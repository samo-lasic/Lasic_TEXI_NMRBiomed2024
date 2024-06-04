function seq = make_seq_FEXI_gc(gc_array, delta_c, delta_pi, delta_max, delta, bf, b_array, tm_array);
gmr = 26.75e7;

if delta > delta_max
    seq = [];
    display('ERROR: delta > delta_max')
    return
end

tm_array = [tm_array(1) tm_array];
tmp_array = tm_array + 2*delta_c/3;
%t_pi = 3/2*delta_pi;
qc_array = gmr*gc_array*delta_c;

Nc = length(gc_array);
Ntm = length(tm_array);
Nb = length(b_array);

seq.n = Nc*Ntm*Nb;
seq.nc = Nc;
seq.ntm = Ntm-1;
seq.nb = Nb;

% filter
gf = sqrt(bf./(gmr^2*delta.^2*(2/3*delta+delta_pi)));
% detect
g_array = sqrt(b_array./(gmr^2*delta.^2*(2/3*delta+delta_pi)));

seq.delta_c = delta_c;
seq.delta_pi = delta_pi;
seq.delta_max = delta_max;
seq.tm_array = tm_array(2:end);
seq.delta = delta;
seq.gf = gf;
seq.bf = bf;
seq.g2_array = g_array;
seq.gc_array = gc_array;
seq.qc_array = qc_array;

c = 0;
for nc = 1:Nc
    for ntm = 1:Ntm
        for nb = 1:Nb
            c = c+1;
            seq.gc(c) = gc_array(nc);
            seq.qc(c) = qc_array(nc);
            seq.gc_ind(c) = nc;

            % mixing
            %bm_array = qc^2*tmp_array;
            seq.bm(c) = seq.qc(c)^2*tmp_array(ntm);

            if ntm == 1
                seq.b1(c) = 0;
                seq.g1(c) = 0;
            else
                seq.b1(c) = bf;
                seq.g1(c) = gf;
            end

            seq.b2(c) = b_array(nb);
            seq.g2(c) = g_array(nb);
            %seq.bm(c) = bm_array(ntm);
            seq.tm(c) = tm_array(ntm);
            %seq.delta(c) = delta;

            seq.b_ind(c) = nb;

        end
    end
end

seq.b = seq.b1+seq.b2+seq.bm;

utm = unique(seq.tm);
for c = 1:length(utm)
seq.tm_ind(seq.tm == utm(c)) = c;
end

seq.s_ind = ones(1,seq.n);
u = unique([seq.b1'*1e-9 seq.tm'*1e3],'rows');
for c = 1:size(u,1)
    ind = find(prod(u(c,:) == [seq.b1'*1e-9 seq.tm'*1e3],2));
    seq.s_ind(ind) = c;
end

end


function seq = make_seq_Vw_gc(gc_array, delta_c, delta_pi, b_array, tm_array, Vw, inlude_DDE)
% seq.b = seq.b1 + seq.bm + seq.b2
% inlude DDE at short tm to minimize exchange sensitivity

gmr = 26.75e7;
tmp_array = tm_array + 2*delta_c/3;
t_pi = 3/2*delta_pi;
qc_array = gmr*gc_array*delta_c;

Nc = length(qc_array);
Ntm = length(tm_array);
Nb = length(b_array);

for nc = 1:Nc
    qc = qc_array(nc);

    for nb = 1:Nb
        b = b_array(nb);
        %delta^2 + delta*t_pi = 3*(b - qc^2*tmp)./(b*Vw - 2*qc^2/delta_c);
        delta2 = 3*(b - qc.^2.*tmp_array)./(b*Vw - 2*qc.^2/delta_c);
        delta = (- t_pi + sqrt(t_pi^2 + 4*delta2) )/2;

        delta_array(nc,:,nb) = delta;
        g2 = 1/(2*gmr)^2*(b*Vw-2*qc.^2/delta_c)./delta;
        g = sqrt(g2);
        g_array(nc,:,nb) = g;

        b1_array(nc,:,nb) = .5*(b-qc^2.*tmp_array);
        bm_array(nc,:,nb) = qc^2.*tmp_array;
    end

end


seq.n = Nc*Nb*(Ntm + inlude_DDE);
seq.nc = Nc;
seq.ntm = Ntm + inlude_DDE;
seq.nb = Nb;

seq.inludes_DDE = inlude_DDE;
seq.delta_c = delta_c;
seq.Vw = Vw;
seq.delta_pi = delta_pi;

c = 0;
for nc = 1:Nc
    for ntm = 1:Ntm

        if inlude_DDE == 1 & ntm == 1
            M = 2;
        else
            M = 1;
        end

        for m = M:-1:1
            for nb = 1:Nb
                c = c+1;
                seq.gc(c) = gc_array(nc);
                seq.gc_ind(c) = nc;
                seq.qc(c) = qc_array(nc);
                seq.b(c) = b_array(nb);
                seq.b_ind(c) = nb;

                seq.isDDE(c) = m - 1;
                seq.b1(c) = b1_array(nc,ntm,nb) * m;
                seq.g1(c) = g_array(nc,ntm,nb) * sqrt(m);
                seq.bm(c) = bm_array(nc,ntm,nb);
                seq.b2(c) = seq.b1(c) * (2-m);
                seq.g2(c) = seq.g1(c) * (2-m);
                %seq.g(c) = g_array(nc,ntm,nb);
                seq.tm(c) = tm_array(ntm);
                seq.tm_ind(c) = ntm;
                seq.delta(c) = delta_array(nc,ntm,nb);
            end
        end
    end
end

seq.gc_array = gc_array;
seq.qc_array = qc_array;
seq.b_array = b_array;
seq.delta_array = delta_array;
% seq.g_array = g_array;

seq.tm_array = tm_array;
if inlude_DDE
    seq.tm_array = [seq.tm_array(1) seq.tm_array];
    seq.tm_ind(~seq.isDDE) = seq.tm_ind(~seq.isDDE) + 1;
end

seq.delta_max = max(seq.delta);

end
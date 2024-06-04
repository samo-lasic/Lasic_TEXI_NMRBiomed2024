function seq = make_seq_b_gc(gc_array, delta_c, delta_pi, delta_max, b_array, tm_array, delta_array, inlude_DDE)
%delta_max determines transverse time (not necessary equal to max of delta_array)
% seq.b = seq.b1 + seq.bm + seq.b2

gmr = 26.75e7;
tmp_array = tm_array + 2*delta_c/3;
t_pi = 3/2*delta_pi;
qc_array = gmr*gc_array*delta_c;

Nc = length(qc_array);
Ntm = length(tm_array);
Nb = length(b_array);
Ndelta = length(delta_array);

for nc = 1:Nc
    qc = qc_array(nc);

    for ntm = 1:Ntm
        tmp = tmp_array(ntm);

        for nb = 1:Nb
            b = b_array(nb);

            g = sqrt( 3/4/gmr^2 * (b - qc^2.*tmp)./delta_array.^2./(delta_array + 3/2*delta_pi) );
            g_array(nc,ntm,nb,:) = g;

            %b1 = (gmr*g).^2.*delta_array.^2.*(2/3*delta_array+delta_pi);
            %b1_array(nc,nb,ntm) = b1(1);

            b1_array(nc,ntm,nb) = .5*(b-qc^2.*tmp);
            bm_array(nc,ntm,nb) = qc^2.*tmp;

        end
    end
end

%bOut = 2*b1_array+b2_array;

seq.n = Nc*Nb*Ndelta*(Ntm + inlude_DDE);
seq.nc = Nc;
seq.ntm = Ntm + inlude_DDE;
seq.nb = Nb;
seq.ndelta = Ndelta;

seq.inludes_DDE = inlude_DDE;

seq.delta_c = delta_c;
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
                for ndelta = 1:Ndelta
                    c = c+1;

                    seq.gc(c) = gc_array(nc);
                    seq.gc_ind(c) = nc;
                    seq.qc(c) = qc_array(nc);
                    seq.b(c) = b_array(nb);
                    seq.b_ind(c) = nb;

                    %                     seq.b1(c) = b1_array(nc,ntm,nb);
                    %                     seq.bm(c) = bm_array(nc,ntm,nb);
                    %                     seq.b2(c) = seq.b1(c);

                    seq.isDDE(c) = m - 1;
                    seq.b1(c) = b1_array(nc,ntm,nb) * m;
                    seq.g1(c) = g_array(nc,ntm,nb,ndelta) * sqrt(m);
                    seq.bm(c) = bm_array(nc,ntm,nb);
                    seq.b2(c) = seq.b1(c) * (2-m);
                    seq.g2(c) = seq.g1(c) * (2-m);

                    % seq.g(c) = g_array(nc,ntm,nb,ndelta);

                    seq.tm(c) = tm_array(ntm);
                    seq.tm_ind(c) = ntm;

                    seq.delta(c) = delta_array(ndelta);
                end
            end
        end
    end
end

seq.gc_array = gc_array;
seq.qc_array = qc_array;
seq.b_array = b_array;

seq.tm_array = tm_array;
seq.g1_array = g_array;
seq.g2_array = g_array;

if inlude_DDE
    seq.tm_array = [seq.tm_array(1) seq.tm_array];
    seq.tm_ind(~seq.isDDE) = seq.tm_ind(~seq.isDDE) + 1;

    tmp_g1_array = seq.g1_array(:,1,:,:) * sqrt(2);
    tmp_g2_array = seq.g2_array(:,1,:,:) * 0;

    tmp_g1_array(:,2,:,:) = seq.g1_array(:,1,:,:);
    tmp_g2_array(:,2,:,:) = seq.g2_array(:,1,:,:);

    tmp_g1_array(:,3:Ntm+1,:,:) = seq.g1_array(:,2:Ntm,:,:);
    tmp_g2_array(:,3:Ntm+1,:,:) = seq.g2_array(:,2:Ntm,:,:);

    seq.g1_array = tmp_g1_array;
    seq.g2_array = tmp_g2_array;
end

seq.delta_array = delta_array;

seq.delta_max = delta_max;
end
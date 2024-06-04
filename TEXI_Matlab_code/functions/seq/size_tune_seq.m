function seq_path = size_tune_seq(R_tune, fn_seq2tune, tol, constraints, tune_step, show_ADC, show_refine_steps, show_tunable_in_full_range)
% size tune sequence
% includes constraints

load(fn_seq2tune)
load(strrep(fn_seq2tune,'.mat','_DvsR.mat'))

display(sprintf('delta max = %g ms', seq.delta_max*1e3));

[~, nR_tune] = min(abs(DvsR.R_array - R_tune));
R_tune = DvsR.R_array(nR_tune);

ADC = DvsR.ADC(:,nR_tune);
minADC = min(ADC);
maxADC = max(ADC);
display(sprintf('ADC range = %g - %g',min(ADC), max(ADC)))

sub_seq.n = seq.n;
sub_seq.nc = seq.nc;
sub_seq.nb = seq.nb;
sub_seq.ntm = seq.ntm;
sub_seq.ndelta = seq.ndelta;

% used to check constraints in find_tunable
sub_seq.g1 = reshape(seq.g1, seq.ndelta, seq.nb, seq.ntm, seq.nc);
sub_seq.g2 = reshape(seq.g2, seq.ndelta, seq.nb, seq.ntm, seq.nc);

ADC = reshape(ADC, seq.ndelta, seq.nb, seq.ntm, seq.nc) / DvsR.D0;
% ADC = ADC / DvsR.D0;

minADC = minADC / DvsR.D0;
maxADC = maxADC / DvsR.D0;

% interpolate delta
Nint = 1e5;
delta_int = interp1(linspace(0,1,seq.ndelta),seq.delta_array,linspace(0,1,Nint));

% figure(1),clf
% hold on
% plot(linspace(0,1,Ndelta),seq.delta_array,'o')
% plot(linspace(0,1,Nint),delta_int,'.')

% g_int = interp1(linspace(0,1,seq.n),seq.g,linspace(0,1,Nint));
% figure(1),clf
% hold on
% plot(linspace(0,1,seq.n),seq.g,'o')
% plot(linspace(0,1,Nint),g_int,'.')


%     tuneADC_array = linspace(0,0.3,1000); % clinical
%     n_tunable = find_tunable(seq, Nint, ADC, tuneADC_array, tol);
%     figure(1),clf
%     plot(tuneADC_array,n_tunable,'.')

% find common ADC

if (show_tunable_in_full_range) % just to check full range
    % tuneADC_array = linspace(0.031,0.033,100); % preclinical R = 2.5 um

    tuneADC_array = linspace(minADC,maxADC,100);

    n_tunable = find_tunable(sub_seq, Nint, ADC, tuneADC_array, tol, constraints);
    figure(1),clf
    plot(tuneADC_array,n_tunable,'.')

    [max_n_tunable, ind] = max(n_tunable);
    tuneADC = tuneADC_array(ind(1));
end


%         plot(ADC(:))
N = 20;
tuneADC_array = linspace(minADC,maxADC,N);
range_ADC_array = 1;
while range_ADC_array > tune_step * (maxADC-minADC)
    [tuneADC_array range_ADC_array, tuneADC] = refine_tuneADC_array(tuneADC_array, N, sub_seq, Nint, ADC, tol, constraints, show_refine_steps);
end

if (show_ADC)
    figure(3),clf
    hold on
    plot(ADC(:))
    % plot([0 seq.n],squeeze(ADC(end,1,1))*[1 1],'r-')
    % plot([0 seq.n],squeeze(ADC(1,end,1))*[1 1],'g-')
    plot([0 seq.n],tuneADC*[1 1],'r--','LineWidth',2)
    %         return
end

delta_tuned = zeros(seq.nc, seq.ntm, seq.nb);
g1_tuned = delta_tuned;
g2_tuned = delta_tuned;
b1_tuned = delta_tuned;
b2_tuned = delta_tuned;
ADCtuned = delta_tuned;
for nc = 1:seq.nc
    for ntm = 1:seq.ntm
        for nb = 1:seq.nb

            g1_int = interp1(linspace(0,1,seq.ndelta),squeeze(seq.g1_array(nc,ntm,nb,:)),...
                linspace(0,1,Nint));

            g2_int = interp1(linspace(0,1,seq.ndelta),squeeze(seq.g2_array(nc,ntm,nb,:)),...
                linspace(0,1,Nint));

            dADC = (ADC/tuneADC-1) * 100;
            dADC_int1 = interp1(linspace(0,1,seq.ndelta),dADC(:,nb, ntm, nc),linspace(0,1,Nint));
            [val, i] = min(abs(dADC_int1));
            if val < tol
                delta_tuned(nc,ntm,nb) = delta_int(i);
   
                g1_tuned(nc,ntm,nb) = g1_int(i);
                g2_tuned(nc,ntm,nb) = g2_int(i);
            end

            ind = find(seq.gc_ind == nc & seq.tm_ind == ntm & seq.b_ind == nb);

            b1_tuned(nc,ntm,nb) = seq.b1(ind(1));
            b2_tuned(nc,ntm,nb) = seq.b2(ind(1));
        end
    end
end

% save sequence
seq_tuned.n = seq.nc * seq.ntm * seq.nb;
seq_tuned.delta_c = seq.delta_c;
seq_tuned.delta_pi = seq.delta_pi;
seq_tuned.nc = seq.nc;
seq_tuned.ntm = seq.ntm;
seq_tuned.nb = seq.nb;

c = 0;
for nc = 1:seq.nc
    for ntm = 1:seq.ntm
        for nb = 1:seq.nb

            c = c+1;
            seq_tuned.gc(c) = seq.gc_array(nc);
            seq_tuned.qc(c) = seq.qc_array(nc);
            seq_tuned.gc_ind(c) = nc;

            seq_tuned.b(c) = seq.b_array(nb);
            seq_tuned.b_ind(c) = nb;

            seq_tuned.tm(c) = seq.tm_array(ntm);
            seq_tuned.tm_ind(c) = ntm;

            seq_tuned.delta(c) = delta_tuned(nc,ntm,nb);
            seq_tuned.g1(c) = g1_tuned(nc,ntm,nb);
            seq_tuned.g2(c) = g2_tuned(nc,ntm,nb);

            seq_tuned.b1(c) = b1_tuned(nc,ntm,nb);
            seq_tuned.b2(c) = b2_tuned(nc,ntm,nb);

        end
    end
end


seq_tuned.delta_max = seq.delta_max;
seq_tuned.b_array = seq.b_array;
seq_tuned.gc_array = seq.gc_array;
seq_tuned.qc_array = seq.qc_array;
seq_tuned.tm_array = seq.tm_array;
seq_tuned.delta_array = delta_tuned;

% clean up - remove non-tunable
ind = find(seq_tuned.delta);
display(sprintf('%d out of %d sequences are tunable',length(ind),seq_tuned.n))

seq_tuned.n = length(ind);
seq_tuned.gc = seq_tuned.gc(ind);
seq_tuned.gc_ind = seq_tuned.gc_ind(ind);
seq_tuned.qc = seq_tuned.qc(ind);
seq_tuned.b = seq_tuned.b(ind);
seq_tuned.b_ind = seq_tuned.b_ind(ind);
seq_tuned.b1 = seq_tuned.b1(ind);
seq_tuned.b2 = seq_tuned.b2(ind);
seq_tuned.tm = seq_tuned.tm(ind);
seq_tuned.tm_ind = seq_tuned.tm_ind(ind);
seq_tuned.bm = seq_tuned.qc.^2.*(seq_tuned.tm+2/3*seq_tuned.delta_c);
seq_tuned.delta = seq_tuned.delta(ind);
seq_tuned.g1 = seq_tuned.g1(ind);
seq_tuned.g2 = seq_tuned.g2(ind);

seq = seq_tuned;
seq.R_tune = R_tune;

display(sprintf('tuned delta = %.2f-%.2f ms, tuned g1 = %.2f-%.2f mT/m, g2 = %.2f-%.2f mT/m',...
    1e3*min(seq.delta),1e3*max(seq.delta), ...
    1e3*min(seq.g1),1e3*max(seq.g1), ...
    1e3*min(seq.g2),1e3*max(seq.g2)))

[folder, name] = fileparts(fn_seq2tune);
name = strrep(name,'b_', sprintf('ST_d%d_',R_tune*2e6));

seq_path = fullfile(folder, name);
save(seq_path,'seq')
display(sprintf('saved: %s', seq_path))

end


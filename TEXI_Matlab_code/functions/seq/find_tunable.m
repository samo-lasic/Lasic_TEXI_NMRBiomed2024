function n_tunable = find_tunable(seq, Nint, ADC, tuneADC_array, tol, constraints)
% find how many sequences are tunable given tuneADC_array thresholds
% includes gradient constraint

parfor n = 1:length(tuneADC_array)
%     for n = 1:length(tuneADC_array)

    tuneADC = tuneADC_array(n);
    dADC = (ADC/tuneADC-1) * 100; %relative difference from target

    n_tunable(n) = 0;

    delta_tuned = zeros(seq.nc, seq.ntm, seq.nb);
    g_tuned = delta_tuned;
    ADCtuned = delta_tuned;
    for nc = 1:seq.nc
        for ntm = 1:seq.ntm
            for nb = 1:seq.nb
                dADC_int1 = interp1(linspace(0,1,seq.ndelta),dADC(:,nb, ntm, nc),linspace(0,1,Nint));
                [val, i] = min(abs(dADC_int1));
                if val < tol & max(seq.g1(:,nb, ntm, nc)) < constraints.g_lim & max(seq.g2(:,nb, ntm, nc)) < constraints.g_lim
                    n_tunable(n) = n_tunable(n) + 1;
                end
            end
        end
    end
end
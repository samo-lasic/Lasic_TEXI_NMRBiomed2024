function aux = Ning_Gamma_aux(seq, isFEXI, dt)
% exchange sensitivity for Ning fit with Gamma approximation

for c = 1:seq.n
    if isFEXI
        [g, q] = def_FEXI(seq.delta, seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);
    else
        [g, q] = def_sequence(seq.delta(c), seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);
    end

    if (0)
        figure(1),clf
        subplot(2,1,1), plot(g)
        subplot(2,1,2), plot(q)
        q(end)
    end

    q2 = q.^2;

    l = length(q2);
    q4 = xcorr(q2,q2) * dt;
    q4 = [0; q4((l+1):end)];

    t = [0:dt:(l-1)*dt]';
    Gamma(c) = 2/seq.b(c)^2 * sum(q4 .* t) * dt;

end

aux.Gamma = Gamma;
aux.tm_ind = seq.tm_ind;
aux.b = seq.b;

end
function Vw = seq2Vw(seq, isFEXI, dt)
% exchange sensitivity 
gmr = 26.75e7;


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

    b = sum(q.^2) * dt;
    Vw(c) = gmr^2 * sum(g.^2) * dt / b;

end


end
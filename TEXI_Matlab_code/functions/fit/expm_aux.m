function aux = expm_aux(seq, isFEXI, dt)
% dephasing waveforms for expm

% find max length
if isFEXI
    g = def_FEXI(seq.delta, seq.delta_max, 0, 0, max(seq.tm), 0, seq.delta_c, seq.delta_pi, dt);
else
    g = def_sequence(max(seq.delta), seq.delta_max, 0, 0, max(seq.tm), 0, seq.delta_c, seq.delta_pi, dt);
end
Lg = length(g);
q2_array = zeros(Lg,seq.n);

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
    l = length(q);
    q2_array(1:l,c) = q.^2;
end

% t = [0:dt:(Lg-1)*dt]';
% aux.t = t;
aux.n = seq.n;
aux.dt = dt;
aux.q2 = q2_array;

aux.tm_ind = seq.tm_ind;
if isfield(seq,'s_ind')
    aux.s_ind = seq.s_ind;
else
    aux.s_ind = seq.tm_ind; % for VT the tm_ind is used instead of s_ind
end
aux.b = seq.b;

end
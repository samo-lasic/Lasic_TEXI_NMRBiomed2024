function aux = Ning_h_aux(seq, isFEXI, dt)
% exchange sensitivity for Ning fit

% find max length
if isFEXI
    g = def_FEXI(seq.delta, seq.delta_max, 0, 0, max(seq.tm), 0, seq.delta_c, seq.delta_pi, dt);
else
    g = def_sequence(max(seq.delta), seq.delta_max, 0, 0, max(seq.tm), 0, seq.delta_c, seq.delta_pi, dt);
end
Lg = length(g);
q4_array = zeros(Lg,seq.n);

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
    q4_array(1:l,c) = q4;
end

t = [0:dt:(Lg-1)*dt]';
aux.t = t;
aux.dt = dt;
aux.q4 = q4_array;

aux.tm_ind = seq.tm_ind;
if isfield(seq,'s_ind')
    aux.s_ind = seq.s_ind;
else
    aux.s_ind = seq.tm_ind; % for VT the tm_ind is used instead of s_ind
end
aux.b = seq.b;

end
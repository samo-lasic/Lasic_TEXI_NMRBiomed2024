function wfm = make_wfm(seq_path, dt);

load(seq_path)
[~, seq_name] = fileparts(seq_path);
isFEXI = contains(seq_name,{'FEXI'});
display(sprintf('isFEXI = %d \n %s', isFEXI, seq_path))

wfm = {};
wfm.n = seq.n;
wfm.dt = dt;

for c = 1:seq.n
    if isFEXI
        [g, q] = def_FEXI(seq.delta, seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);
    else
        [g, q] = def_sequence(seq.delta(c), seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);
    end

    wfm.seq(c).g = g;
    wfm.seq(c).q = q;

    if (0)
        figure(1),clf
        subplot(2,1,1), plot(g)
        subplot(2,1,2), plot(q)
    end

end



function save_seq_with_const_b(fn_ref_seq, fn_out_seq, constraints, inlude_DDE, Ndelta, dt, seqName_out, show_b)
% fn_ref_seq is a references sequence, e.g. with constant VT from which we take some parameters

% read from Vw = const.
load(fn_ref_seq)

delta_min = constraints.delta_min; %min(seq.delta);
delta_max = max(seq.delta);
display(sprintf('delta = %g-%g ms',delta_min*1e3,delta_max*1e3))

delta_array = linspace(delta_min,delta_max, Ndelta);

seq = make_seq_b_gc(seq.gc_array, seq.delta_c, seq.delta_pi, delta_max, seq.b_array, seq.tm_array, delta_array, inlude_DDE);

for c = 1:seq.n
    [g, q] = def_sequence(seq.delta(c), seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);
    b(c) = sum(q.^2)*dt;

    %     figure(1),clf
    %     subplot(2,1,1), plot(g)
    %     subplot(2,1,2), plot(q)
    %     pause(.1)   

    m0Test(c) = abs(q(end))/max(abs(q)); % just checking m0 = 0 requirement
end

if max(m0Test) > 1e-3
    display('problems with m0!')
end


if (show_b)
    figure(1),clf
    hold on
    plot(seq.b,'o')
    plot(b,'.')
end

mkdir(fileparts(fn_out_seq));

save(fn_out_seq,'seq')
display(sprintf('saved: %s', fn_out_seq))


display(sprintf('max g1/g2 = %g/%g', max(seq.g1 * 1000), max(seq.g2 * 1000)))

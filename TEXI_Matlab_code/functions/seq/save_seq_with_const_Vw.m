function seq = save_seq_with_const_Vw(seq_folder, seq_name, gc_array, delta_c, delta_pi, b_array, tm_array, Vw, inlude_DDE, do_save)
% max gradient determines bmin (via delta) and Vw

seq = make_seq_Vw_gc(gc_array, delta_c, delta_pi, b_array, tm_array, Vw, inlude_DDE);

% b-value = red-blue, crusher = dark-bright

str{1} = sprintf('gap180 = %.0f ms, gc = %.0f - %.0f mT/m, delta_c = %.2f ms, Vw = %g s^{-2}',...
    seq.delta_pi*1e3, min(seq.gc_array)*1e3, max(seq.gc_array)*1e3, ...
    seq.delta_c*1e3,seq.Vw);

str{2} = sprintf('delta = %.2f-%.2f ms, g = %.2f-%.2f mT/m, b = %.2f-%.2f s/mm^2', ...
    min(seq.delta)*1e3, max(seq.delta)*1e3,...
    min(seq.g1)*1e3, max(seq.g1)*1e3, min(seq.b_array)*1e-6,max(seq.b_array)*1e-6);

figure(2),clf
subplot(2,1,1)
%col = jet(seq.nb);
col = [linspace(1,0,seq.nb)' 0*linspace(1,0,seq.nb)' linspace(0,1,seq.nb)'];
hold on
for nc = 1:seq.nc
    if seq.nc > 1 
        brightness = (seq.gc_array(nc)-min(seq.gc))/(max(seq.gc)-min(seq.gc));
    else
        brightness = 1;
    end
    for nb = 1:seq.nb
        c = brightness*col(nb,:);

        ind = find(seq.gc_ind == nc & seq.b_ind == nb);
        delta_array = seq.delta(ind) * 1e3;
        tm_array = seq.tm(ind) * 1e3;

        plot(tm_array, delta_array,'.-','color',c, 'MarkerSize',18)
    end
end
ylabel('delta [ms]')
set(gca,'FontSize',14)
box off

title(str{1},'FontSize',10)


subplot(2,1,2)
hold on
for nc = 1:seq.nc
    if seq.nc > 1
        brightness = (seq.gc_array(nc)-min(seq.gc))/(max(seq.gc)-min(seq.gc));
    else
        brightness = 1;
    end

    for nb = 1:seq.nb
        c = brightness*col(nb,:);
        %g_array = squeeze(seq.g_array(nc,:,nb))'*1e3;
        ind = find(seq.gc_ind == nc & seq.b_ind == nb);
        g_array = max([seq.g1(ind);seq.g2(ind)]) * 1e3;
        tm_array = seq.tm(ind) * 1e3;
        plot(tm_array, g_array,'.-','color',c, 'MarkerSize',18)
    end
end

ylabel('g_{max} [mT/m]')
xlabel('tm [ms]')
set(gca,'FontSize',14)
box off

title(str{2},'FontSize',10)


display_seq_par_range(seq_name, seq, 0);

if do_save
    mkdir(seq_folder)
    seq_path = fullfile(seq_folder, seq_name);
    save(seq_path,'seq')
    display(sprintf('saved: %s', seq_path))
else
    display(sprintf('not saving: %s', seq_name))
end

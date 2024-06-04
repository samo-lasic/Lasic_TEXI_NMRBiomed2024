function fh = show_seq_par(n_fig, fn, isFEXI, seq, opt)

if nargin < 5
    opt.lw = 4;
    opt.ms = 60;
    opt.fs = 24;
    opt.show_all = 0;
    opt.show_labels = 1;
end

if isFEXI
    display_seq_par_range(fn, seq, isFEXI);

    fh = figure(n_fig);
    clf;
    hold on

    b = seq.b1+seq.b2;   
    bf_ind = find(seq.b1 == 0);
    b_ind = find(seq.b1 > 0);
    
    
    plot(1e3*seq.tm(bf_ind),1e-6*b(bf_ind),'.r','LineWidth',opt.lw,'MarkerSize',1.5*opt.ms);
    plot(1e3*seq.tm(b_ind),1e-6*b(b_ind),'.k','LineWidth',opt.lw,'MarkerSize',opt.ms);

    set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)

    limX = 1e3*[min(seq.tm) max(seq.tm)];
    limY = 1e-6*[0 max([seq.b1 + seq.b2])];

    xlim(limX + .1*diff(limX)*[-1 1])
    ylim(limY + .1*diff(limY)*[-1 1])

    if opt.show_labels
        ylabel('b [ms/mm^2]')
        xlabel('tm [ms]')
        title(fn,'interpreter','none','FontSize',opt.fs*.9)
    end

    pbaspect([1 1 1])


else
    display_seq_par_range(fn, seq, isFEXI);

    ind = seq.b2 == 0;

    fh = figure(n_fig);
    clf;
    hold on
    plot(1e3*seq.tm(ind), 1e-6*seq.b(ind),'.r','LineWidth',opt.lw,'MarkerSize',1.5*opt.ms);
    plot(1e3*seq.tm(~ind), 1e-6*seq.b(~ind),'.k','LineWidth',opt.lw,'MarkerSize',opt.ms);
    

    if opt.show_all
        if isfield(seq,'b1')
            plot(1e3*seq.tm, 1e-6*seq.b1,'xr','LineWidth',opt.lw,'MarkerSize',opt.ms/2);
        end

        if isfield(seq,'bm')
            plot(1e3*seq.tm, 1e-6*seq.bm,'^g','LineWidth',opt.lw,'MarkerSize',opt.ms/3);
        end

        if isfield(seq,'b2')
            plot(1e3*seq.tm, 1e-6*seq.b2,'ob','LineWidth',opt.lw,'MarkerSize',opt.ms/3);
        end
    end
    pbaspect([1 1 1])

    limX = 1e3*[min(seq.tm) max(seq.tm)];
    limY = 1e-6*[0 max([seq.b])];

    xlim(limX + .1*diff(limX)*[-1 1])
    ylim(limY + .1*diff(limY)*[-1 1])


    set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)
    if opt.show_labels
        ylabel('b [ms/mm^2]')
        xlabel('t_m [ms]')
        title(fn,'interpreter','none','FontSize',opt.fs*.9)
    end

end
end
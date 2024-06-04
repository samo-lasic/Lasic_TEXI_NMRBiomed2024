function fh = plot_k_list(opt)
opt.nR = 1;


fh = figure(opt.fig_n);
clf
set(gcf, 'InvertHardCopy', 'off')
set(fh,'Color','white', 'Units', 'inches', ...
    'PaperPosition',[0 0 opt.size_mm / 25.4 ], 'PaperPositionMode', 'auto');
screen_size = get(0,'screensize');
fh.Position(3:4) = opt.display_ratio * opt.size_mm;

hold on

cnt = 0;
for n = 1:numel(opt.res_paths)

    % load fit results
    load(opt.res_paths{n},'res')
    isFEXI = contains(opt.seq_names{n},{'FEXI'});

    nc_ind = [];
    if isempty(opt.gc.show)
        nc_ind = 1:res.in.Nc;
    else
        if opt.gc.show(1)
            nc_ind = [nc_ind 1];
        end
        if opt.gc.show(2)
            nc_ind = [nc_ind res.in.Nc];
        end
    end

%     if ~isempty(opt.gc.colormap)
%         opt.gc.col = eval(sprintf('%s(%d)',opt.gc.colormap, res.in.Nc));
%     end

    Xplot_array = res.in.k_array;
    opt.XLIM = 1.02*max(Xplot_array)*[0 1];


    switch opt.plot_case
        case 1 % k vs k loop:crushers

            opt.plot_name = 'k';
            opt.Xlabel = 'true k [1/s]';
            opt.Ylabel = 'estimated k [1/s]';
            Yplot_array = squeeze(res.fitP_array(:,:,:,3));

            opt.YLIM = opt.XLIM;
            %opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];


        case 2 % V/MD2 or sigma vs k loop:crushers

            opt.plot_name = 'V_k';
            opt.Xlabel = 'true k [1/s]';

            if isFEXI
                Yplot_array = squeeze(res.fitP_array(:,:,:,2));
                opt.Ylabel = 'sigma';
            else
                Yplot_array = squeeze(res.fitP_array(:,:,:,2) ./ res.fitP_array(:,:,:,1).^2);
                opt.Ylabel = 'V / MD^2';
            end

            opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];


        case 3 % MD or ADC0 vs k loop:crushers

            opt.plot_name = 'MD_k';
            opt.Xlabel = 'true k [1/s]';

            if isFEXI
                Yplot_array = 1e9*squeeze(res.fitP_array(:,:,:,1));
                opt.Ylabel = 'ADC_0 [mm^2/ms]';
            else
                Yplot_array = 1e9*squeeze(res.fitP_array(:,:,:,1));
                opt.Ylabel = 'MD [mm^2/ms]';
            end

            opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];


        case 4 % residual sum of squares
            if ~isfield(res,'sum_square_res')
                display('fit_model_for_nonFEXI_data not present in res')
                fh = [];
                return
            end

            opt.plot_name = 'rss';
            opt.Xlabel = 'true k [1/s]';

            Yplot_array = 1e6*res.sum_square_res;
            opt.Ylabel = 'RSS [x10^{-6}]';
            %opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];
            opt.YLIM = [0 2*mean(Yplot_array(:))] .* [.95 1.05];

    end

   

    for nc = nc_ind
        Yplot = squeeze(Yplot_array(nc, :, opt.nR));

        colormap = eval(sprintf('%s(%d)',opt.gc.colormap, res.in.Nc));
        col = colormap(nc,:);

        if isfield(opt,'substrate')
            line_style = opt.substrate.line_style(n);
            lw = opt.substrate.lw(n);
            marker = opt.substrate.marker(n);
            ms = opt.substrate.ms(n);

            if isfield(opt.substrate,'colormap') % override color
                colormap = eval(sprintf('%s(%d)',opt.substrate.colormap, numel(opt.res_paths)));
                col = colormap(n,:);
            end
        else
            line_style = opt.gc.line_style(nc);
            lw = opt.gc.lw(nc);
            marker = opt.gc.marker(nc);
            ms = opt.gc.ms(nc);

            
        end

        ph = plot(Xplot_array, Yplot, 'LineStyle',line_style,...
            'Marker', marker, 'MarkerSize',ms,...
            'LineWidth',lw); %,'Color',opt.col(n,:))


        if isfield(opt,'col')
            ph.Color = opt.col(n,:);
        else
            ph.Color = col;
        end

        cnt = cnt+1;

        if opt.show_empty_str
            legend_str{cnt} = ' ';
        else
            legend_str{cnt} = [opt.legend_str{n} '-c' num2str(round(res.in.qc_array(nc)))]; %num2str(nc>1)];
        end
    end
end

if  opt.plot_case == 1
    plot(opt.XLIM,opt.XLIM,'--k')
    if isempty(find(opt.XLIM - opt.YLIM))
        axis equal
    end
end



if opt.show_legend
    lh = legend(legend_str,'interpreter','none','location','northwest','FontSize',opt.legend_fs, ...
        'Box','on', 'Color','white','EdgeColor','white');
end

set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out',...
    'TickLength',[.02 .02],'FontSize',opt.fs)

xlim(opt.XLIM)
ylim(opt.YLIM)
pbaspect([1 1 1])
    
if opt.plot_case == 1
    grid on
    grid minor
    a = .03;
    set(gca,'XTick',get(gca,'YTick'), 'GridAlpha', 2*a,'MinorGridAlpha',a,'MinorGridLineStyle','-');
%     axis equal
end


if opt.show_label
    xlabel(opt.Xlabel)
    ylabel(opt.Ylabel)
end


end





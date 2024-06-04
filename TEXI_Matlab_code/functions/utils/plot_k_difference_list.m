function fh = plot_k_difference_list(opt)
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

    % nc_ind = [];
    % if opt.gc.show(1)
    %     nc_ind = [nc_ind 1];
    % end
    % if opt.gc.show(2)
    %     nc_ind = [nc_ind res.in.Nc];
    % end

    nc_ind = [1 res.in.Nc];

    Xplot_array = res.in.k_array;
    opt.XLIM = 1.1*max(Xplot_array)*[0 1];


    switch opt.plot_case
        case 1 % k vs k loop:crushers

            opt.plot_name = 'dif_k';
            opt.Xlabel = 'true k [1/s]';
            opt.Ylabel = 'dif. k fit (w/o crush) [1/s]';
            Yplot_array = squeeze(res.fitP_array(:,:,:,3));

            opt.YLIM = opt.XLIM;
            %opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];


        case 2 % V/MD2 or sigma vs k loop:crushers

            opt.plot_name = 'dif_V_k';
            opt.Xlabel = 'true k [1/s]';

            if isFEXI
                Yplot_array = squeeze(res.fitP_array(:,:,:,2));
                opt.Ylabel = 'dif. sigma (w/o crush) ';
            else
                Yplot_array = squeeze(res.fitP_array(:,:,:,2) ./ res.fitP_array(:,:,:,1).^2);
                opt.Ylabel = 'dif. V / MD^2 (w/o crush) ';
            end

            opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];


        case 3 % MD or ADC0 vs k loop:crushers

            opt.plot_name = 'dif_MD_k';
            opt.Xlabel = 'true k [1/s]';

            if isFEXI
                Yplot_array = 1e9*squeeze(res.fitP_array(:,:,:,1));
                opt.Ylabel = 'dif. ADC_0 (w/o crush) [mm^2/ms]';
            else
                Yplot_array = 1e9*squeeze(res.fitP_array(:,:,:,1));
                opt.Ylabel = 'dif. MD (w/o crush) [mm^2/ms]';
            end

            opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];


        case 4 % residual sum of squares
            if ~isfield(res,'sum_square_res')
                display('fit_model_for_nonFEXI_data not present in res')
                fh = [];
                return
            end

            opt.plot_name = 'dif_rss';
            opt.Xlabel = 'true k [1/s]';

            Yplot_array = 1e6*res.sum_square_res;
            opt.Ylabel = 'dif. RSS (w/o crush) [x10^{-6}]';
            %opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];
            opt.YLIM = [0 2*mean(Yplot_array(:))] .* [.95 1.05];

    end


%     if isfield(opt,'substrate')
%         line_style = opt.substrate.line_style(n);
%         lw = opt.substrate.lw(n);
%     else
%         line_style = opt.gc.line_style(nc);
%         lw = opt.gc.lw(nc);
%     end


    Yplot = squeeze(Yplot_array(nc_ind(2), :, opt.nR)) - squeeze(Yplot_array(nc_ind(1), :, opt.nR));


          if isfield(opt,'substrate')
            line_style = opt.substrate.line_style(n);
            lw = opt.substrate.lw(n);
            marker = opt.substrate.marker(n);
            ms = opt.substrate.ms(n);
        else
            line_style = opt.gc.line_style(2);
            lw = opt.gc.lw(2);
            marker = opt.gc.marker(2);
            ms = opt.gc.ms(2);
        end

        ph = plot(Xplot_array, Yplot, 'LineStyle',line_style,...
            'Marker', marker, 'MarkerSize',ms,...
            'LineWidth',lw); %,'Color',opt.col(n,:))

%     ph = plot(Xplot_array, Yplot, 'LineStyle',line_style,...
%         'Marker', opt.gc.marker{2}, 'MarkerSize',opt.gc.marker_size(2),...
%         'LineWidth',lw); %,'Color',opt.col(n,:))


    %     plot(Xplot_array, Yplot, 'LineStyle',opt.gc.line_style{2},...
    %         'Marker', opt.gc.marker{2}, 'MarkerSize',opt.gc.marker_size(2),...
    %         'LineWidth',opt.gc.lw,'Color',opt.col(n,:))

    if ~isempty(opt.gc.colormap)
        opt.gc.col = eval(sprintf('%s(%d)',opt.gc.colormap, res.in.Nc));
    end


    if isempty(opt.gc.colormap)
        ph.Color = opt.col(n,:);
    else
        ph.Color = opt.gc.col(2,:);
    end


    cnt = cnt+1;

    legend_str{cnt} = opt.legend_str{n}; % '-c' num2str(round(res.in.qc_array(nc)))]; %num2str(nc>1)];

end

% if  opt.plot_case == 1
%     plot(opt.XLIM,opt.XLIM,'--k')
%     if isempty(find(opt.XLIM - opt.YLIM))
%         axis equal
%     end
% end
xlim(opt.XLIM)
% ylim(opt.YLIM)
pbaspect([1 1 1])



lh = legend(legend_str,'interpreter','none','Box','off','location','best','FontSize',opt.legend_fs);

% title_str{1} = sprintf('%s: qc = %.0f m^{-1}, R = %.1f - %.1f \\mum', ...
%     strrep(extractAfter(res.in.fn,'_'),'_','-'), qc, min_R_array*1e6, max_R_array*1e6);


set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out',...
    'TickLength',[.02 .02],'FontSize',opt.fs)

%     grid on
%     grid minor

xlabel(opt.Xlabel)
ylabel(opt.Ylabel)


end


% function fh = plot_loop_crusher(fig_n_start, nc_ind, Xplot_array, Yplot_array, res, opt)
%
% for nc = nc_ind
%     %gc = res.in.gc_array(nc);
%     qc = res.in.qc_array(nc);
%     min_R_array = min(res.in.R_array);
%     max_R_array = max(res.in.R_array);
%
%     % fitP_names = {'MD','VD','k'};
%     % fitP_name = fitP_names{opt.n_fit_par};
%     %k_fit_array = squeeze(res.fitP_array(nc,:,:,3));
%
%     col = [linspace(1,0,res.in.NR)' 0*linspace(0,1,res.in.NR)' linspace(0,1,res.in.NR)'];
%
%     fh(nc) = figure(fig_n_start+nc-1);
%     clf
%     hold on
%     for nR = 1:res.in.NR
%         Yplot = squeeze(Yplot_array(nc, :, nR));
%         plot(Xplot_array,Yplot,'.-','MarkerSize',16,'LineWidth',1,'Color',col(nR,:))
%     end
%
%     if  opt.plot_case == 1
%         plot(opt.XLIM,opt.XLIM,'--k')
%         if isempty(find(opt.XLIM - opt.YLIM))
%             axis equal
%         end
%     end
%     xlim(opt.XLIM)
%     ylim(opt.YLIM)
%     pbaspect([1 1 1])
%
%     title_str{1} = sprintf('%s: qc = %.0f m^{-1}, R = %.1f - %.1f \\mum', ...
%         strrep(extractAfter(res.in.fn,'_'),'_','-'), qc, min_R_array*1e6, max_R_array*1e6);
%
%     if isfield(res.in,'fg') & isfield(res.in,'Dg')
%         title_str{2} = sprintf('fg = %.2f, Dg = %g m^2/s', ...
%             res.in.fg, res.in.Dg);
%     end
%
%     set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out','TickLength',[.02 .02],...
%         'FontSize',opt.fs)
%
% %     grid on
% %     grid minor
%
%     xlabel(opt.Xlabel)
%     ylabel(opt.Ylabel)
%
%     %title(title_str,'interpreter','none')
%
%     th = title(title_str);
%     if isfield(opt, 'title_relative_fs')
%         th.FontSize = th.FontSize * opt.title_relative_fs;
%     end
%     if isfield(opt, 'title_color')
%         th.Color = opt.title_color;
%     end
%
%     if opt.save_fig
%
%         figName = [res.in.fn '_' opt.plot_name '_gc' num2str(nc) '.png' ];
%
%         if isfield(opt,'save_path')
%             path = opt.save_path;
%         else
%             path = fullfile(fileparts(pwd), 'results/fig');
%         end
%
%         mkdir(path)
%         path = fullfile(path, figName);
%         display(sprintf('saving ... %s',path))
%         print(fh(nc),path,'-dpng','-r300');
%     end
% end
%
%
% end




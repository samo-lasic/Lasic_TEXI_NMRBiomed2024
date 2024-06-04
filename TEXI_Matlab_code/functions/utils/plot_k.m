function fh = plot_k(fig_n_start, res, opt)


% Ning_h [MD V k S0]
% FEXI11 [ADC0 sigma AXR S0]

% res.fitP_array [Nc, Nk, NR, Npar]

if nargin < 3
    opt.lw = 3;
    opt.fs = 14;
    opt.ms = 12;
    opt.save_fig = 0;
    opt.plot_case = 1;
    opt.plot_extreme_crushers = 0;
end

isFEXI = contains(res.in.fn,{'FEXI'});

if opt.plot_extreme_crushers
    nc_ind = [1 res.in.Nc];
else
    nc_ind = 1:res.in.Nc;
end

Xplot_array = res.in.k_array;
opt.XLIM = 1.1*max(Xplot_array)*[0 1];

switch opt.plot_case
    case 1 % k vs k loop:crushers

        opt.plot_name = 'k';
        opt.Xlabel = 'k [1/s]';
        opt.Ylabel = 'k fit [1/s]';
        Yplot_array = squeeze(res.fitP_array(:,:,:,3));

        % opt.YLIM = opt.XLIM;
        opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];

        fh = plot_loop_crusher(fig_n_start, nc_ind, Xplot_array, Yplot_array, res, opt);

    case 2 % V/MD2 or sigma vs k loop:crushers

        opt.plot_name = 'V_k';
        opt.Xlabel = 'k [1/s]';

        if isFEXI & str2num(extractAfter(res.in.fn,'_fit')) == 1
            Yplot_array = squeeze(res.fitP_array(:,:,:,2));
            opt.Ylabel = 'sigma';
        else
            Yplot_array = squeeze(res.fitP_array(:,:,:,2) ./ res.fitP_array(:,:,:,1).^2);
            opt.Ylabel = 'V / MD^2';
        end

        opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];

        fh = plot_loop_crusher(fig_n_start, nc_ind, Xplot_array, Yplot_array, res, opt);

    case 3 % MD or ADC0 vs k loop:crushers

        opt.plot_name = 'MD_k';
        opt.Xlabel = 'k [1/s]';

        if isFEXI & str2num(extractAfter(res.in.fn,'_fit')) == 1
            Yplot_array = 1e9*squeeze(res.fitP_array(:,:,:,1));
            opt.Ylabel = 'ADC_0 [mm^2/ms]';
        else
            Yplot_array = 1e9*squeeze(res.fitP_array(:,:,:,1));
            opt.Ylabel = 'MD [mm^2/ms]';
        end

        opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];

        fh = plot_loop_crusher(fig_n_start, nc_ind, Xplot_array, Yplot_array, res, opt);

    case 4 % residual sum of squares
        if ~isfield(res,'sum_square_res')
            display('fit_model_for_nonFEXI_data not present in res')
            fh = [];
            return
        end

        opt.plot_name = 'rss';
        opt.Xlabel = 'k [1/s]';

        Yplot_array = 1e6*res.sum_square_res;
        opt.Ylabel = 'RSS [x10^{-6}]';
        %opt.YLIM = [min(Yplot_array(:)) max(Yplot_array(:))] .* [.95 1.05];
        opt.YLIM = [0 2*mean(Yplot_array(:))] .* [.95 1.05];

        fh = plot_loop_crusher(fig_n_start, nc_ind, Xplot_array, Yplot_array, res, opt);



end

end


function fh = plot_loop_crusher(fig_n_start, nc_ind, Xplot_array, Yplot_array, res, opt)

for nc = nc_ind
    %gc = res.in.gc_array(nc);
    qc = res.in.qc_array(nc);
    min_R_array = min(res.in.R_array);
    max_R_array = max(res.in.R_array);

    % fitP_names = {'MD','VD','k'};
    % fitP_name = fitP_names{opt.n_fit_par};
    %k_fit_array = squeeze(res.fitP_array(nc,:,:,3));

    col = [linspace(1,0,res.in.NR)' 0*linspace(0,1,res.in.NR)' linspace(0,1,res.in.NR)'];

    
    if fig_n_start > 0
        fh(nc) = figure(fig_n_start+nc-1);
    else
        fh(nc) = figure(length(findobj('type','figure')) + 1);
    end
    fh(nc).Color = 'white';

    clf
    hold on
    for nR = 1:res.in.NR
        Yplot = squeeze(Yplot_array(nc, :, nR));
        plot(Xplot_array,Yplot,'.-','MarkerSize',16,'LineWidth',1,'Color',col(nR,:))
    end

    if  opt.plot_case == 1
        plot(opt.XLIM,opt.XLIM,'--k')
        if isempty(find(opt.XLIM - opt.YLIM))
            axis equal
        end
    end
    xlim(opt.XLIM)
    ylim(opt.YLIM)
    pbaspect([1 1 1])

%     title_str{1} = sprintf('%s: \n qc = %.0f m^{-1}, R = %.1f - %.1f \\mum', ...
%         strrep(extractAfter(res.in.fn,'_'),'_','-'), qc, min_R_array*1e6, max_R_array*1e6);

    title_str{1} = sprintf('%s: \n qc = %.0f m^{-1}, R = %.1f - %.1f \\mum', ...
        strrep(res.in.fn,'_','-'), qc, min_R_array*1e6, max_R_array*1e6);

    if isfield(res.in,'fg') & isfield(res.in,'Dg')
        title_str{2} = sprintf('fg = %.2f, Dg = %g m^2/s', ...
            res.in.fg, res.in.Dg);
    end

    set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out','TickLength',[.02 .02],...
        'FontSize',opt.fs)

%     grid on
%     grid minor

    xlabel(opt.Xlabel)
    ylabel(opt.Ylabel)

    %title(title_str,'interpreter','none')

    th = title(title_str);
    if isfield(opt, 'title_relative_fs')
        th.FontSize = th.FontSize * opt.title_relative_fs;
    end
    if isfield(opt, 'title_color')
        th.Color = opt.title_color;
    end 

    if opt.save_fig

        figName = [res.in.fn '_' opt.plot_name '_gc' num2str(nc) '.png' ];

        if isfield(opt,'save_path')
            path = opt.save_path;
        else
            path = fullfile(fileparts(pwd), 'results/fig');
        end

        mkdir(path)
        path = fullfile(path, figName);
        display(sprintf('saving ... %s',path))
        print(fh(nc),path,'-dpng','-r300');
    end
end


end




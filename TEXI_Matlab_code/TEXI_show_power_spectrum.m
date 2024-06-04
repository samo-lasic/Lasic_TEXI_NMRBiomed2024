% help script to plot power spectra

clear all
close all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

seq_folder = fullfile(root_dir, 'seq_TEXI');

fn_append = '';  % ignore if empty
protocols = {};
protocols(end+1).seq_name = 'FEXI_preclinX_dz1000_tm400_bf2000_b1300';
protocols(end).tm_ind = [1 3]; % (FEXI) if empty do all
protocols(end).filter_names = {'bmin'}; % ignore if empty
protocols(end+1).seq_name = 'VT_preclinX_dz1000_tm250_b2600';
protocols(end).tm_ind = [2 6]; % (TEXI) if empty do all
protocols(end).filter_names = {}; % ignore if empty


plot_wfm = 1;
plot_spectra = 1;
plot_attenuation = 1;

opt.save_fig = 0;

% tm_ind = [2 6]; % if empty do all
gc_ind = [1 4]; % if empty do all
b_ind = []; % if empty do max b

show_min_exchange_weighting = 0; % show sequences with minimum exchange weighting, i.e. only filter block in FEXI, only one transverse encoding block in TEXI

dt = 1e-6;
N0pad = 1e7;
fThresh = 1000;


opt.fig_n = 0; % ofset figure number
opt.fig_width_mm = 30;
opt.fig_aspect = .09;
opt.display_ratio = .15;

opt.fLIM = [0 100];
opt.lw1 = 4;
opt.lw2 = 1.5;
opt.ms = 10;
opt.fs = 20; %20;
opt.title_fs_ratio = .5;
% opt.show_all = 0;
opt.show_labels = 0;
opt.show_yticklabels = 1;

all_seq_files = dir(seq_folder);


for n_protocol = 1:numel(protocols)
    seq_name = protocols(n_protocol).seq_name;
    tm_ind = protocols(n_protocol).tm_ind; % if empty do all
    filter_names = protocols(n_protocol).filter_names; % ignore if empty

    ind = ones(1, numel(all_seq_files));

    if ~isempty(fn_append)
        ind = ind & contains({all_seq_files.name}, fn_append);
    end

    % optionally add seq name filter here
    if ~isempty(seq_name)
        ind = ind & contains({all_seq_files.name}, seq_name);

        if ~isempty(filter_names)
            ind = ind & ~contains({all_seq_files.name}, filter_names);
        end
    end

    ind = find(ind & ~contains({all_seq_files.name},{'_DvsR','b_'}));
    seq_file = all_seq_files(ind);
    if numel(seq_file) > 1
        display('multiple sequence files')
        continue
    end


    [~,seq_name,~] = fileparts(seq_file.name);

    load(fullfile(seq_file.folder,seq_file.name))

    isFEXI = contains(seq_file.name,{'FEXI_'});

    % select sequence parameters
    if isempty(tm_ind)
        sub_ind = ones(1,length(seq.tm_ind));
    else
        sub_ind = ismember(seq.tm_ind, tm_ind);
    end

    if ~isempty(gc_ind)
        sub_ind = sub_ind & ismember(seq.gc_ind, gc_ind);
    end

    if isempty(b_ind)
        sub_ind = sub_ind & seq.b_ind == max(seq.b_ind);
    else
        sub_ind = sub_ind & ismember(seq.b_ind, b_ind);
    end

    if ~show_min_exchange_weighting
        if isFEXI
            sub_ind = sub_ind & seq.b1 > 0;
        else
            sub_ind = sub_ind & ~seq.isDDE;
        end
    end

    sub_ind = find(sub_ind);
    % sub_ind = find(ismember(seq.tm_ind, tm_ind)  & ismember(seq.gc_ind, gc_ind) & seq.b_ind == max(seq.b_ind));

    sub_seq = get_sub_seq_from_ind(seq, sub_ind);


    Tenc_max = 2*(max(sub_seq.delta) + sub_seq.delta_pi + sub_seq.delta_max + sub_seq.delta_c) + max(sub_seq.tm);

    Nmax = round(Tenc_max/dt);
    NFT = Nmax + N0pad;
    f = linspace(-1,1,NFT)'/2/dt;

    flim = fThresh*[0 1];
    ind_PS = find(f >= flim(1) & f<= flim(2));
    fPS = f(ind_PS);


    % get waveforms & spectra
    qs = zeros(Nmax,sub_seq.n);
    s = [];
    for c = 1:sub_seq.n
        if isFEXI
            [g, q] = def_FEXI(sub_seq.delta, sub_seq.delta_max, sub_seq.g1(c), sub_seq.g2(c), sub_seq.tm(c), sub_seq.gc(c), sub_seq.delta_c, sub_seq.delta_pi, dt);
        else
            [g, q] = def_sequence(sub_seq.delta(c), sub_seq.delta_max, sub_seq.g1(c), sub_seq.g2(c), sub_seq.tm(c), sub_seq.gc(c), sub_seq.delta_c, sub_seq.delta_pi, dt);
        end
        s(:,c) = power_spectrum(q, dt, NFT, ind_PS);
        qs(1:length(q),c) = q;
    end

    protocols(n_protocol).seq_file = seq_file;
    protocols(n_protocol).qs = qs;
    protocols(n_protocol).maxq = max(qs(:));
    protocols(n_protocol).s = s;
    protocols(n_protocol).maxs = max(s(:));

    protocols(n_protocol).Nmax = Nmax;
    protocols(n_protocol).Tenc_max = Tenc_max;
    protocols(n_protocol).fPS = fPS;
    protocols(n_protocol).sub_seq = sub_seq;
    %     protocols(n_protocol).Nmax = Nmax;

end
maxq = max([protocols.maxq]);
maxs = max([protocols.maxs]);

opt.size_mm = opt.fig_width_mm * [1 opt.fig_aspect * sub_seq.n];

for n_protocol = 1:numel(protocols)
    p = protocols(n_protocol);

    if plot_wfm
        % plot waveforms
        % color_order = turbo(sub_seq.n);
        color_order = 0*[0 0 1; 1 0 0; 0 0 1; 1 0 0;];

        t = 1e3*linspace(0,p.Tenc_max,p.Nmax);

        for c = 1:p.sub_seq.n

            fh = figure;
            clf
            set(gcf, 'InvertHardCopy', 'off')
            set(fh,'Color','white', 'Units', 'inches', ...
                'PaperPosition',[0 0 opt.size_mm / 25.4 ], 'PaperPositionMode', 'auto');
            screen_size = get(0,'screensize');
            fh.Position(3:4) = opt.display_ratio * opt.size_mm;

            Y = p.qs(:,c)/maxq;

            X = t;
            X = interp1(1:length(X),X, linspace(1,length(X),2e7));
            Y = interp1(1:length(Y),Y, linspace(1,length(Y),2e7));

            plot(X,Y,'.k','LineWidth',opt.lw1,'MarkerSize',opt.ms,'color',color_order(c,:));
            xlim([0 1e3*p.Tenc_max])

            ylim([0 1])
            if ~opt.show_yticklabels
                yticklabels([])
            end

            title_str{c} = sprintf('tm = %d ms qc = %g m^{-1}', round(p.sub_seq.tm(c) * 1e3), p.sub_seq.qc(c));

            if opt.show_labels
                th = title(title_str{c});
                th.FontSize = th.FontSize * opt.title_fs_ratio;
                if c == 1
                    ylabel('q(t) [a.u.]')
                elseif c == sub_seq.n
                    xlabel('time [ms]')
                end
            end
            set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)


            if opt.save_fig

                fig_name = strrep(p.seq_file.name,'.mat','_wfm');
                fig_name = [fig_name num2str(c)];
                fig_path = fullfile(p.seq_file.folder, 'fig');
                mkdir(fig_path)
                fig_path = fullfile(fig_path, fig_name);

                display(sprintf('saving ... %s',fig_path))
                print(fh, fig_path,'-dpng','-r300');
            end

        end
    end

    % plot spectra
    if plot_spectra
        % color_order = turbo(sub_seq.n);
        color_order = 0*[0 0 1; 1 0 0; 0 0 1; 1 0 0;];

        for c = 1:p.sub_seq.n

            fh = figure;
            clf
            set(gcf, 'InvertHardCopy', 'off')
            set(fh,'Color','white', 'Units', 'inches', ...
                'PaperPosition',[0 0 opt.size_mm / 25.4 ], 'PaperPositionMode', 'auto');
            screen_size = get(0,'screensize');
            fh.Position(3:4) = opt.display_ratio * opt.size_mm;

            X = p.fPS;
            Y = p.s(:,c)/maxs;

            X = interp1(1:length(X),X, linspace(1,length(X),2e7));
            Y = interp1(1:length(Y),Y, linspace(1,length(Y),2e7));

            plot(X,Y,'.k','LineWidth',opt.lw1,'MarkerSize',opt.ms,'color',color_order(c,:));
            if ~isempty(opt.fLIM)
                xlim(opt.fLIM)
            end
            ylim([0 1])
            if ~opt.show_yticklabels
                yticklabels([])
            end

            if opt.show_labels
                %title_str = sprintf('tm = %d ms qc = %g m^{-1}', round(sub_seq.tm(c) * 1e3), sub_seq.qc(c));
                th = title(title_str{c});
                th.FontSize = th.FontSize * opt.title_fs_ratio;

                if c == 1
                    ylabel('power [a.u.]')
                elseif c == sub_seq.n
                    xlabel('frequency [Hz]')
                end
            end
            set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)


            if opt.save_fig

                fig_name = strrep(p.seq_file.name,'.mat','_spec');
                fig_name = [fig_name num2str(c)];
                fig_path = fullfile(p.seq_file.folder, 'fig');
                mkdir(fig_path)
                fig_path = fullfile(fig_path, fig_name);

                display(sprintf('saving ... %s',fig_path))
                print(fh, fig_path,'-dpng','-r300');
            end

        end

    end

    % illustrate tuning with cumulative attenuation
    if plot_attenuation

        fig_width_mm = opt.fig_width_mm * 1.1;
        fig_aspect = 0.9;
        display_ratio = opt.display_ratio;
        size_mm = fig_width_mm * [1 fig_aspect];

        fs = opt.fs;
        lw1 = 4;
        ms = 10;
        lw2 = 2;
        grid_alpha = 0.2;

        title_fs_ratio = 0.8;


        D0 = 2e-9;
        R_array = 0.5*[1:2:11] * 1e-6;
        color_order = redblue(length(R_array));
        color_order = fliplr(color_order);

        %     Dw = fDwCylinder((2*pi*fPS).^2,R_array,D0,0,50)'; % cylindrical restriction
        Dw = fDwSphere((2*pi*p.fPS).^2,R_array,D0,0,50); % spherical restriction

        %Dw is proportional to G * R^4/D0 * w^2
        Dnorm = GeometryFactor(3,100) * R_array.^4 / D0;
        Dw = Dw ./ repmat(Dnorm,length(p.fPS),1);

        clear betas
        for n = 1:p.sub_seq.n
            b = sum(s(:,n));
            betas(n,:,:) = cumsum(p.s(:,n) .* Dw);
        end
        max_betas = max(betas(:));
        YLIM = [0 1.05 * max_betas];


        for c = 1:p.sub_seq.n

            beta = squeeze(betas(c,:,:));

            if (1)
                XLIM = [10 1000]; %[];
                xscale = 'log';
            else
                XLIM = [0 500]; %[];
                xscale = 'lin';
            end

            fh = figure;
            clf
            set(gcf, 'InvertHardCopy', 'off')
            set(fh,'Color','white', 'Units', 'inches', ...
                'PaperPosition',[0 0 size_mm / 25.4 ], 'PaperPositionMode', 'auto');
            screen_size = get(0,'screensize');
            fh.Position(3:4) = display_ratio * size_mm;

            X = p.fPS;
            Y = beta;

            plot(X,Y,'.','LineWidth',lw1,'MarkerSize',ms)

            if ~isempty(YLIM)
                ylim(YLIM)
            end
            if ~isempty(XLIM)
                xlim(XLIM)
            end
            
            grid on
            % grid minor
            set(gca, 'GridColor',[0 .0 .0],'GridLineStyle','-','GridAlpha',grid_alpha);
            set(gca,'XScale', xscale,'YScale', 'lin', 'LineWidth',lw2,'Box','on','TickDir','out','TickLength',[.02 .02],'FontSize',fs)
            % box off

            if mod(c,2) == 0
                set(gca, 'XDir', 'reverse','YAxisLocation', 'right','GridColor',[0 .0 .0],'GridLineStyle','-','GridAlpha',grid_alpha);
            end

            if ~opt.show_yticklabels
                yticklabels([])                
            end


            if opt.show_labels
                th = title(title_str{1});
                th.FontSize = th.FontSize * title_fs_ratio;
            end
            colororder(color_order);

            if opt.save_fig

                fig_name = strrep(p.seq_file.name,'.mat','_beta');
                fig_name = [fig_name num2str(c)]
                fig_path = fullfile(p.seq_file.folder, 'fig');
                mkdir(fig_path)
                fig_path = fullfile(fig_path, fig_name);

                display(sprintf('saving ... %s',fig_path))
                print(fh, fig_path,'-dpng','-r300');
            end
        end

    end

end



% ------------------------------------ FUNCTIONS ---------------------------
function G = GeometryFactor(geometry, Nterms)
% geometry = 1,2,3 (plane R = d/2, cylinder R, sphere R)

if geometry == 1
    % d = 2*R
    % B = 8*d^2./(2*k-1).^4/pi^4;
    % a = (2*k-1).^2*pi^2/d^2;

    k = 1:Nterms;

    K2 = (2*k-1).^2;
    C1 = 8./K2/pi^2;
    C2 = K2.^2*pi^4;

    G = sum(C1./C2);
    G = G * 2^4; % adjusted for R = d/2

elseif geometry == 2
    % B=2*(R./Kc).^2./(Kc.^2-1);
    % a=(Kc/R).^2;

    Kc = BesselKernelsCylinder(Nterms);

    Kc4 = Kc.^4;
    C =2./(Kc.^2-1);
    G = sum(C./Kc4);


elseif geometry == 3
    % B=2*(R./Ks).^2./(Ks.^2-2);
    % a=(Ks/R).^2;

    Ks = BesselKernelsSphere(Nterms);

    Ks4 = Ks.^4;
    C =2./(Ks.^2-2);
    G = sum(C./Ks4);

else
    G = 0;
end

end

function f=BesselKernelsSphere(n)
for l=1:n
    f(l)=fzero(inline('x*besselj(1/2,x)-2*besselj(3/2,x)'),2+(l-1)*pi);
end
end

function g=BesselKernelsCylinder(n)
for l=1:n
    g(l)=fzero(inline('besselj(0,x)-besselj(1,x)/x'),2+(l-1)*pi);
end
end


function help_plot(X,Y,lw1,lw2,ms,grid_alpha, xscale,fs, XLIM, YLIM, direction)
plot(X,Y,'.','LineWidth',lw1,'MarkerSize',ms)
ylim(YLIM)
if ~isempty(XLIM)
    xlim(XLIM)
end
if ~opt.show_yticklabels
    yticklabels([])
end
grid on
% grid minor
set(gca, 'GridColor',[0 .0 .0],'GridLineStyle','-','GridAlpha',grid_alpha);
set(gca,'XScale', xscale,'YScale', 'lin', 'LineWidth',lw2,'Box','on','TickDir','out','TickLength',[.02 .02],'FontSize',fs)
% box off

if direction < 0
    set(gca, 'XDir', 'reverse','YAxisLocation', 'right','GridColor',[0 .0 .0],'GridLineStyle','-','GridAlpha',grid_alpha);
end

end

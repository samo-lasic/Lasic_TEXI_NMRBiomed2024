% show sequence parameters

clear all
close all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

seq_folder = fullfile(root_dir, 'seq');

fn_append = '_preclin7';
select_seq_names = {'FEXI_','VT_'};
%select_seq_names = {'FEXI_','VT_preclin7_dz1000_tm250_b2600_reduceBrange'};

seq_folder = fullfile(root_dir, 'seq_TEXI'); 
select_seq_names = {'FEXI_','VT_preclin7_dz1000_tm250_b2600.'};
%select_seq_names = {'VT_preclin7_dz1000_tm250_b2600','VT_preclin7_dz1000_tm250_b2600_bmin1450','VT_preclin7_dz1000_tm250_b2600_bmin1834'};

% IVIM
% fn_append = '_IVIM';
% select_seq_names = {'FEXI_'};


opt.save_fig = 0;

% -----------------------------------------------------
%show_parameters= [gc   b   g   delta   tm  Gamma   Vw];
show_parameters = [1    1   1   1       1   1       1];

max_gc = 60*1e-3; %[] % if empty show max

delta_range = [2 10]*1e-3; %[] % if empty set auto
tm_range = [0 400]*1e-3; %[] % if empty set auto
Gamma_range = [0 200]*1e-3; %[] % if empty set auto
one_over_sqrtVw_range = .8*1e-3; %[3.7 6.2]*1e-3; %[] % if empty set auto


dt = 1e-4;

nFig = 0; % ofset figure number


files = dir(seq_folder);
ind = contains({files.name}, fn_append); %sprintf('%s_%s_', seqName_in, fn_append));

% optionally add seq name filter here
if ~isempty(select_seq_names)
    ind = ind & contains({files.name},select_seq_names);
end

ind = find(ind & ~contains({files.name},{'_DvsR','b_','_b0'}));
files = files(ind);


for nSeq = 1:numel(files)
    seq_file = files(nSeq);

    [~,seq_name,~] = fileparts(seq_file.name);

    isFEXI = contains(seq_name,{'FEXI_'});

    XLIM = [];
    load(fullfile(seq_file.folder,seq_file.name))

    if isempty(max_gc)
        [~, seq] = get_sub_seq(seq, [1 seq.nc]);
    else
        [~,ind] = min(abs(seq.gc-max_gc)); 
         ind = seq.gc_ind(ind);
         [~, seq] = get_sub_seq(seq, [1 ind(1)]);
    end
 

    opt.lw = 4;
    opt.fs = 24;
    opt.ms = 60;
    opt.show_all = 0;
    opt.show_labels = 0;

    % show only key parameters
    fh = show_seq_par(nSeq, seq_name, isFEXI, seq, opt);

    fh.Color = [1 1 1];

    plot_ind = find(show_parameters);
    n_plots = numel(plot_ind);


    if opt.save_fig

        fig_name = strrep(seq_file.name,'.mat','_main_par');
        fig_path = fullfile(seq_file.folder, 'fig'); 
        mkdir(fig_path)
        fig_path = fullfile(fig_path, fig_name);

        display(sprintf('saving ... %s',fig_path))
        print(fh, fig_path,'-dpng','-r300');
    end


    if (1) % show parameters per sequence
        %opt.ms = 16;
        opt.fs = 16;
        opt.ms1 = 10;
        opt.ms2 = 6;
        opt.ms3 = 16; 
        
        opt.lw1 = 1;
        opt.lw2 = 2;
        opt.show_all = 1;
        opt.show_labels = 0;

        fh = figure(numel(files) + nSeq + nFig);
        clf
        fh.Color = [1 1 1];
        fh.Position(4) = fh.Position(4) * n_plots/3;

        for n_plot = 1:n_plots
            subplot(n_plots,1,n_plot)
            if opt.show_labels & n_plot == 1
                title(seq_name,'interpreter','none','FontSize',opt.fs*.9)
            end

            switch plot_ind(n_plot) % %show_parameters= [gc b    g   delta   tm  Gamma];

                case 1 % gc
                    plot(seq.gc*1e3,'.k','LineWidth',opt.lw1,'MarkerSize',opt.ms3);
                    if opt.show_labels
                        ylabel('g_c [mT/m]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)

                case 2 % b
                    hold on
                    plot(seq.b*1e-6,'.k','LineWidth',opt.lw1,'MarkerSize',opt.ms3);
                    if opt.show_all
                        if isfield(seq,'b1')
                            plot(seq.b1*1e-6,'.r','LineWidth',opt.lw1,'MarkerSize',opt.ms1);
                        end
                        if isfield(seq,'b2')
                            plot(seq.b2*1e-6,'ob','LineWidth',opt.lw1,'MarkerSize',opt.ms2);
                        end
                    end

                    if opt.show_labels
                        ylabel('b [s/mm^2]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)

                case 3 % g
                    hold on

                    if opt.show_all
                        if isfield(seq,'g1')
                            plot(seq.g1*1e3,'.r','LineWidth',opt.lw1,'MarkerSize',opt.ms1);
                        end
                        if isfield(seq,'g2')
                            plot(seq.g2*1e3,'ob','LineWidth',opt.lw1,'MarkerSize',opt.ms2);
                        end
                    else
                        plot(max([seq.g1;seq.g2])*1e3,'.-k','LineWidth',opt.lw1,'MarkerSize',opt.ms3);
                    end

                    if opt.show_labels
                        ylabel('g [mT/m]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)

                case 4 % delta
                    if numel(seq.delta) == 1 
                        delta = ones(1,seq.n) * seq.delta;
                    else
                        delta = seq.delta;
                    end

                    plot(delta*1e3,'.k','LineWidth',opt.lw1,'MarkerSize',opt.ms3);

                    if ~isempty(delta_range)
                        ylim(delta_range*1e3)
                    end

                    if opt.show_labels
                        ylabel('delta [ms]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)

                case 5 % tm
                    plot(seq.tm*1e3,'.k', 'LineWidth',opt.lw1,'MarkerSize',opt.ms3);
                    
                     if ~isempty(tm_range)
                        ylim(tm_range*1e3)
                    end
                    
                    if opt.show_labels
                        ylabel('t_m [ms]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)

                case 6 % Gamma
                    Gamma = Ning_Gamma(seq, isFEXI, dt);
                    plot(Gamma*1e3,'.k', 'LineWidth',opt.lw1,'MarkerSize',opt.ms3);

                     if ~isempty(Gamma_range)
                        ylim(Gamma_range*1e3)
                     end

                    if opt.show_labels
                        ylabel('Gamma [ms]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)
                case 7 % Vw
                    Vw = seq2Vw(seq, isFEXI, dt);
                    y = round(100 * 1./sqrt(Vw)*1e3)/100;
                    plot(y,'.k', 'LineWidth',opt.lw1,'MarkerSize',opt.ms3);

                     if ~isempty(one_over_sqrtVw_range)
                        YLIM = mean(y) + [-1/2 1/2] * one_over_sqrtVw_range*1e3;
                        ylim(YLIM)
                     end


                    if opt.show_labels
                        ylabel('1/sqrt(Vw) [ms]')
                    end
                    set(gca,'YScale', 'lin','LineWidth',opt.lw2,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)
            end
            if isempty(XLIM) 
                XLIM = [1 seq.n];
            end

            xlim(XLIM)

        end

        if opt.show_labels
            %            ylabel('t_m [ms]')
            xlabel('sequence')
        end

    end

    if opt.save_fig

        fig_name = strrep(seq_file.name,'.mat','_par');
        fig_path = fullfile(seq_file.folder, 'fig'); 
        mkdir(fig_path)
        fig_path = fullfile(fig_path, fig_name);

        display(sprintf('saving ... %s',fig_path))
        print(fh, fig_path,'-dpng','-r300');
    end



end




% inspect raw signals, fitted signals and fit parameters for selected substrate/size, crusher and exchange rate  

clear all
close all

restoredefaultpath

addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

root_path = fullfile(root_dir, 'data_TEXI', 'results','sim07');
seq_path = fullfile(root_dir, 'seq_TEXI');


% --------

results_filter = {}; % can be empty for no filter
% results_filter = {'_b2600.mat'}

results_filter_inlude = 1; % 1/0 include/exclude



sim_folder_name = 'sim07';
substrate_name = 'signals_spheres_d5_regular_rho50';
seq_names = {'FEXI_preclin7_dz500_tm400_bf2000_b1300', 'VT_preclin7_dz1000_tm250_b2600'};
seq_names = seq_names(2);
results_filter = {'fit3_b2600'}



plot_fit_parameters = 1;
plot_signals = 1; % show fit signals
plot_raw_signals = 0; % for checking raw signals (script has pause!)

plot_signal_for_d = 5*1e-6; % if only one size is in res.in.R_array it will find that size

plot_signals_for_all_k = 0; % if 0 plot only extreme cases (depends on plot_signals_for_selected_k)
plot_signals_for_selected_k = []; % if empty use all k else select k
plot_signals_for_all_gc = 2; % if 0 plot only extreme cases


opt.save_fig = 0;

opt.plot_case = 1;
% 1 - k vs k loop:crushers, 2 - V or sigma vs k loop:crushers,
% 3 - MD or ADC0 vs k loop:crushers, 4 - residual sum of squares vs k loop:crushers
opt.plot_extreme_crushers = ~plot_signals_for_all_gc; % if 1 plot only min and max crushers

opt.title_relative_fs = 0.7; % comment this for default

opt.lw = 3;
opt.fs = 14;
opt.ms = 12;

opt.showGrid = 1;

substrate_name = strrep(substrate_name,'signals', sim_folder_name);


for n_fn = 1:numel(seq_names)
    seq_name = seq_names{n_fn};

    % search candidate sequences
    d_seq = dir(seq_path);
    d_seq = d_seq(contains({d_seq.name}, seq_name) & ~contains({d_seq.name}, 'DvsR'));
    seq_name_array = {d_seq.name};
    seq_name_array = strrep(seq_name_array,'.mat','');

    sim_data_names = strcat(seq_name_array,'_',substrate_name);

    % find matching sim_data_names in results folder (containing part of sequence name)
    d_sim = dir(root_path);
    res_files = d_sim(contains({d_sim.name}, sim_data_names) & ~contains({d_sim.name}, 'some string'));

    % filter names
    if ~isempty(results_filter)
        if results_filter_inlude
            res_files = res_files(contains({res_files.name}, results_filter));
        else
            res_files = res_files(~contains({res_files.name}, results_filter));
        end
    end


    colors = turbo(numel(res_files))*.7;
    for n_res = 1:numel(res_files)

        fn = fullfile(res_files(n_res).folder,res_files(n_res).name);


        % load fit results
        try
            load(fn,'res')
            display(sprintf('found: %s', fn))
        catch
            display(sprintf('not found: %s', fn))
            continue
        end


        if (plot_fit_parameters)
            opt.save_path = fullfile(root_path, substrate_name, seq_name);
            opt.title_color = colors(n_res,:);
            fh = plot_k(0, res, opt);
        end

        % --------------------  show signal decays ----------------
        if (plot_signals)


            seq_name = extractBefore(res_files(n_res).name,'_sim');
            isFEXI = contains(seq_name,{'FEXI'});

            seq_file = fullfile(root_dir,'seq',seq_name);
            load(seq_file)


            [~,nR] = min(abs(res.in.R_array - .5*plot_signal_for_d));


            d = round(2e6 * res.in.R_array(nR));

            if plot_signals_for_all_gc
                nc_loop = [1:res.in.Nc];
            else
                nc_loop = [1 res.in.Nc];
            end

            if isempty(plot_signals_for_selected_k)
                sel_k = [1:res.in.Nk];
            else
                sel_k = plot_signals_for_selected_k;
            end
            if plot_signals_for_all_k
                nk_loop = sel_k;
            else
                nk_loop = [sel_k(1) sel_k(end)];
            end


            for nc = nc_loop
                qc = res.in.qc_array(nc);
                for nk = nk_loop
                    k = res.in.k_array(nk);

                    [sub_ind, sub_seq] = get_sub_seq(seq, nc);

                    sFIT = squeeze(res.fit_sig_array(nk,nR,sub_ind));
                    sig = squeeze(res.sig_array(nk,nR,sub_ind));
                    Pout = squeeze(res.fitP_array(nc,nk,nR,:));

                    opt.title_str = sprintf('%s: \n qc = %.1f m^{-1}, d = %d \\mum, k = %.1f s^{-1}', ...
                        strrep(res.in.fn,'_','-'), qc, d, k);

                    opt.title_color = colors(n_res,:);


                    fh = fixi_fig(length(findobj('type','figure')) + 1, isFEXI, sub_seq, sig, sFIT, Pout, opt);

                    if opt.save_fig
                        if isfield(opt,'save_path')
                            path = opt.save_path;
                        else
                            path = root_path;
                        end
                        mkdir(path)
                        figName = [res.in.fn '_' 'sig' '_d' num2str(d) '_k' num2str(k) '_gc' num2str(nc) '.png' ];

                        path = fullfile(path, figName);
                        display(sprintf('saving ... %s',path))
                        print(fh,path,'-dpng','-r300');

                    end

                end
            end
        end


        if (plot_raw_signals)

            seq_name = extractBefore(res_files(n_res).name,'_sim');

            seq_file = fullfile(root_dir,'seq',seq_name);
            load(seq_file)


            [~,nR] = min(abs(res.in.R_array - .5*plot_signal_for_d));
            d = round(2e6 * res.in.R_array(nR));

            figure, clf
            for nc = 1:res.in.Nc
                [sub_ind, sub_seq] = get_sub_seq(seq, nc);
                sig_array = squeeze(res.sig_array(:,nR, sub_ind));
                [~,ind] = sort(repmat(1:sub_seq.nb,1,sub_seq.ntm));

                for nk = 1:res.in.Nk
                    s = sig_array(nk,:);
                    plot(reshape(s(ind), sub_seq.ntm, sub_seq.nb), '-o','linewidth', 2);
                    xlabel('tm - ind')
                    ylabel('signal')
                    title(sprintf('qc = %d, k = %d',res.in.qc_array(nc), res.in.k_array(nk)))
                    pause();
                end
            end


        end

    end

end


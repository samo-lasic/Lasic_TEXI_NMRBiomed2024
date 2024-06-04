
% Plot results for the blood-brain barrier exchange estimated with the original AXR fit and the Ning model.

% combine plots with different substrates and crushers etc.
% assumes each result is for a single substrate (size)
% plotting only min and/or max crusher data
% only k vs k

% see the fitting script for more info regarding subsampling of simulation data

clear all
% close all

restoredefaultpath

addpath(genpath(fullfile(pwd,'functions')));

warning('off','all')


%root_dir = fileparts(fileparts(pwd));
root_dir = pwd;


root_path = fullfile(root_dir, 'res_IVIM');

fig_opt = 0;


sim_cases = 9; %9:12; %1:8;
fit_models = 1:2;

for sim_case = sim_cases
    for fit_model = fit_models

        result_name = sprintf('sim_%d_FEXI_IVIM_fit%d.mat',sim_case,fit_model);
        result_names = {result_name};
        fig_opt = 2;
        

        N_res = numel(result_names);

        opt.save_fig = 1;
        %opt.append_name = ''; % add custom string for saving figures

        opt.plot_case = 1; % which parameter to plot
        % 1 - k vs k loop:crushers, 2 - V or sigma vs k loop:crushers,
        % 3 - MD or ADC0 vs k loop:crushers, 4 - residual sum of squares vs k loop:crushers
        opt.plot_difference = 0;

        opt.show_legend = 0;
        opt.show_empty_str = 0;
        opt.show_label = 0;

        opt.YLIM = [0.1 1]; % comment this for default
        opt.title_relative_fs = 0.9; % comment this for default

        opt.lw = 4;
        opt.fs = 14;

        opt.legend_fs = 11.5;
        opt.size_mm = [1 1] * 6;
        opt.display_ratio = 1;


        % select crusher data to plot
        opt.gc.show = []; %show min [1 0] and/or max [0/1 1] or all crushers []

        opt.gc.line_style = {'-', '-','-','-','-'};
        opt.gc.lw = 2 * ones(1,numel(opt.gc.line_style)); %if ~isfield(opt,'substrate') use gc.lw
        opt.gc.marker = {'o', '.','.','.','.'};
        opt.gc.ms = 20 * ones(1,numel(opt.gc.marker));
        opt.gc.ms(1) = 12;
        opt.gc.colormap = 'copper'; %'viridis';%[];%'copper' %[]; % if [] use opt.col else use colormap for different crushers

        if ~isfield(opt,'append_name')
            opt.append_name = '';
        end

        switch fig_opt % if ~isfield(opt,'substrate') use opt.gc
            case 1

                % select crusher data to plot
                opt.gc.show = []; %show min [1 0] and/or max [0/1 1] or all crushers []

                opt.substrate.line_style = {'-', ':','-', ':'};
                opt.substrate.lw = 2 * ones(1,numel(opt.substrate.line_style));
                opt.substrate.marker = {'o', '.','o', '.'};
                opt.substrate.ms = 20 * ones(1,numel(opt.substrate.marker));
                opt.substrate.ms(1) = 12;
                % opt.substrate.colormap = 'viridis'; %'copper' %if field exists override use this else use gc.colormap, if opt.col exists this will be ignored

            case 2
                % select crusher data to plot
                opt.gc.show = []; %show min [1 0] and/or max [0/1 1] or all crushers []

                ms = 1e-6 + 0; % 1 or 1e-6 to show/hide
                %opt.substrate.line_style = {':', '-'};
                opt.substrate.line_style = {'-', ':','--','-.'};
                opt.substrate.lw = 2 * ones(1,numel(opt.substrate.line_style));
                opt.substrate.marker = {'o', '.','x','*'};
                opt.substrate.ms = 12 * ones(1,numel(opt.substrate.marker)) * ms;
                opt.substrate.ms(2) = 20 * ms;
                %opt.substrate.colormap = 'viridis';%[];%'copper' %[]; % if [] use opt.col else use colormap for different crushers


            case 3
                % select crusher data to plot
                opt.gc.show = []; %show min [1 0] and/or max [0/1 1] or all crushers []

                opt.substrate.line_style = {'-'}; opt.substrate.line_style = repmat(opt.substrate.line_style,1,6);
                opt.substrate.lw = 2 * ones(1,numel(opt.substrate.line_style));
                opt.substrate.marker = {'.'}; opt.substrate.marker = repmat(opt.substrate.marker,1,6);
                opt.substrate.ms = 20 * ones(1,numel(opt.substrate.marker));
                %opt.substrate.colormap = 'viridis';%[];%'copper' %[]; % if [] use opt.col else use colormap for different crushers

            case 4
                % select crusher data to plot
                opt.gc.show = [0 1]; %show min [1 0] and/or max [0/1 1] or all crushers []

                opt.substrate.line_style = {'-'}; opt.substrate.line_style = repmat(opt.substrate.line_style,1,6);
                opt.substrate.lw = 2 * ones(1,numel(opt.substrate.line_style));
                opt.substrate.marker = {'.'}; opt.substrate.marker = repmat(opt.substrate.marker,1,6);
                opt.substrate.ms = 20 * ones(1,numel(opt.substrate.marker));
                opt.substrate.colormap = 'redblue'; %'viridis' 'inferno';%[];%'copper' %[]; % if [] use opt.col else use colormap for different crushers

        end



        % split names, make paths and legend strings
        for n_res = 1:numel(result_names)
            res_name = result_names{n_res};
            res_path = fullfile(root_path, res_name);
            seq_name = 'FEXI'; %extractBefore(res_name,'_sim');
            seq_name_short = extractBefore(res_name,'_dz');
            substrate_name = extractAfter(res_name,'_sim');
            sim_name =  strcat('sim',extractBefore(substrate_name,'_'));
            substrate_name = extractAfter(substrate_name,'_');
            substrate_name = extractBefore(substrate_name, '_fit');
            fit_name = extractAfter(res_name,'_fit');

            opt.sim_names{n_res} = sim_name;
            opt.substrate_names{n_res} = substrate_name;
            opt.seq_names{n_res} = seq_name;
            opt.res_paths{n_res} = res_path;

            opt.legend_str{n_res} = strcat(seq_name_short,'-',substrate_name,'_fit',fit_name);

        end

        substrate_name = find_common_substrings(opt.substrate_names);
        sim_name = find_common_substrings(opt.sim_names);
        seq_name = find_common_substrings(opt.seq_names);

        fig_folder = strcat(sim_name,'_',substrate_name);
        fig_folder = fullfile(root_path,'fig',fig_folder);

        fig_name = strrep(result_name,'.mat','.png');

        opt.title_str = '';

        opt.fig_n = 1;
        fh = plot_k_list(opt);
        if opt.show_label
            title(fig_name,'interpreter','none')
        end

        if opt.save_fig

            mkdir(fig_folder)
            fig_path = fullfile(fig_folder, fig_name);
            display(sprintf('saving ... %s',fig_path))
            print(fh, fig_path,'-dpng','-r300');
        end

        if opt.plot_difference
            opt.fig_n = 2;
            fh = plot_k_difference_list(opt);

            if opt.save_fig

                mkdir(fig_folder)
                fig_path = fullfile(fig_folder, [fig_name '_dif']);
                display(sprintf('saving ... %s',fig_path))
                print(fh, fig_path,'-dpng','-r300');
            end
        end
    end
end

function common_substrings = find_common_substrings(strs)

% Split each string into substrings
all_substrings = cellfun(@(x) split(x, '_'), strs, 'UniformOutput', false);

% Find the common substrings
common_substrings = all_substrings{1};
for i = 2:length(all_substrings)
    common_substrings = intersect(common_substrings, all_substrings{i});
end

% Convert the result to a single string
common_substrings = char(join(common_substrings, '_'));

end
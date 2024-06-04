% make new data and seq with a reduced set of b values for the VT sequence
% all b at short tm and high b at long tm
% output is labeled as '..._reduceBrange'
% output is not labeled depending on bmin - take care not to overwrite!
% see the fitting script for more info regarding subsampling of simulation data

clear all
close all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

simulation_data_folder = fullfile('data_TEXI','simulation');
seq_folder = 'seq_TEXI';


fn_append = 'preclin7';
fn_substrate_tag = '07'; % different substrates for the same sequences / simulation folder

sim_case = 1; % 1-3
switch sim_case
    case 1
        sim_fn_array = {...
            'signals_spheres_d1_regular_rho50',...
            'signals_spheres_d3_regular_rho50',...
            'signals_spheres_d5_regular_rho50',...
            'signals_spheres_d7_regular_rho50',...
            'signals_spheres_d9_regular_rho50',...
            'signals_spheres_d11_regular_rho50'};
        % sim_fn_array = strcat(sim_fn_array,'_rho50_20231024')

    case 2
        sim_fn_array = {...
            'signals_spheres_d4_gamma_dist',...
            'signals_spheres_d8_gamma_dist'};

    case 3
        sim_fn_array = {...
            'signals_cylinders_d1_regular', ...
            'signals_cylinders_d3_regular', ...
            'signals_cylinders_d5_regular', ...
            'signals_cylinders_d7_regular', ...
            'signals_cylinders_d9_regular', ...
            'signals_cylinders_d11_regular'};
        sim_fn_array = strcat(sim_fn_array,'_rho50_20231024')
end



% b_min = 1450;
b_min = 1834;



for n_sim = 1:numel(sim_fn_array)
    sim_fn = sim_fn_array(n_sim);

    if numel(sim_fn) > 1 sim_fn = sim_fn(1); end
    signal_fn = fullfile(root_dir, simulation_data_folder, fn_substrate_tag, char(sim_fn));


    % load simulated signals
    load(signal_fn)
    signal_names = fields(signals);

    % select non FEXI sequences
    isFEXI = contains(signal_names,'FEXI_');
    signal_names = signal_names(~isFEXI);

    clear seq_files

    % get sequence paths from seq folder
    for n = 1:numel(signal_names)
        seq_files(n) = dir(fullfile(root_dir, seq_folder, [signal_names{n} '.mat']));
    end

    % clean up any empty fields
    seq_files = seq_files(find(cellfun(@ischar,{seq_files.name})));
    seq_names = {seq_files.name};

    for n_seq = 1:numel(seq_names)

        seq_name = seq_names{n_seq};

        %load sequence
        file = seq_files(find(contains(seq_names, seq_name)));
        seq_name = file.name;
        load(fullfile(file.folder, seq_name))
        seq_name = extractBefore(seq_name,'.mat');


        % limit low b-value for longer mixing times
        bmin = (.95*b_min*1e6) * ones(size(seq.tm_array));
        bmin(1) = 0;

        % select seqeunces
        sub_ind = seq.tm < 0;
        for i = 1:length(seq.tm_array)
            sub_ind = sub_ind | (seq.tm == seq.tm_array(i) & seq.b > bmin(i));
        end

        sub_ind = find(sub_ind);
        seq = get_sub_seq_from_ind(seq, sub_ind);

%         strAppend = '_reduceBrange'; % this should be fixed so that the name reflects the b-limits
        strAppend = sprintf('_bmin%d',b_min);

        % new sequence and data names
        sub_seq_name = [seq_name strAppend];
        [~,~,ext] = fileparts(signal_fn);
        if isempty(ext)
            new_signal_fn = [signal_fn strAppend];
        else
            new_signal_fn = strrep(signal_fn, ext, strAppend);
        end


        % save new sequence
        save(fullfile(file.folder, sub_seq_name),'seq')

        % make new signals structure
        clear sub_signals

        sub_sig_struct = signals.(seq_name);

        % get all relevant fields
        f = fields(sub_sig_struct);
        f = f(contains(f,'_k'));
        N = numel(f);
        for n = 1:N;
            s = signals.(seq_name).(f{n});
            sub_signals.(sub_seq_name).(f{n}) = s(sub_ind);
        end

    end

    % save new signals
    signals = sub_signals;

    save(new_signal_fn,'signals');

end






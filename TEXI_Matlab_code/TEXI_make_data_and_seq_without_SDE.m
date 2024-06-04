% A sub-sampling script.
% make new data and seq without the SDE in VT and ST sequences 
% append _noSDE
% used in testing


clear all
close all

restoredefaultpath

addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

simulation_data_folder = fullfile('data_TEXI','simulation');
seq_folder = 'seq_TEXI';

fn_append = 'preclin7';
fn_substrate_tag = '07';
sim_fn_array = {'signals_spheres_d3_regular_rho50.mat'};


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

        % select seqeunces with b2 > 0
        sub_ind = find(seq.b2 > 0);
        seq = get_sub_seq_from_ind(seq, sub_ind);

        % save new sequence
        sub_seq_name = [seq_name '_noSDE'];
        save(fullfile(file.folder, sub_seq_name),'seq')

        % make new signals structure
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
    new_signal_fn = strrep(signal_fn,'.mat','_noSDE.mat');
    save(new_signal_fn,'signals');

end






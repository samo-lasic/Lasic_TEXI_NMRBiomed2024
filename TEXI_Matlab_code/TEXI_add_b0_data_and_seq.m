% make new data and seq by additing S = 1 @ b = 0
% appending "_b0"
% for VT data and sequence
% used in testing

clear all
close all

restoredefaultpath

addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

simulation_data_folder = fullfile('data_TEXI','simulation');
seq_folder = 'seq_TEXI';

fn_prepend = '';
fn_append = 'preclin7';
fn_substrate_tag = '07';
sim_fn_array = {'signals_spheres_d3_regular_rho50.mat'};

seq_include = {'VT_'}; % if empty include all else select those containing strings in the list



for n_sim = 1:numel(sim_fn_array)
    sim_fn = sim_fn_array(n_sim);


    if numel(sim_fn) > 1 sim_fn = sim_fn(1); end

    signal_fn = fullfile(root_dir, simulation_data_folder, fn_substrate_tag, char(sim_fn));


    % load simulated signals
    load(signal_fn)
    signal_names = fields(signals);

    % get sequence paths from seq folder
    for n = 1:numel(signal_names)
        seq_files(n) = dir(fullfile(root_dir, seq_folder, [signal_names{n} '.mat']));
    end

    % clean up empty fields
    seq_files = seq_files(find(cellfun(@ischar,{seq_files.name})));
    seq_names = {seq_files.name};

    if ~isempty(seq_include)
        seq_files = seq_files(contains(seq_names,seq_include));
        seq_names = seq_names(contains(seq_names,seq_include));
    end


    for n_seq = 1:numel(seq_names)

        seq_name = seq_names{n_seq};

        %load sequence
        file = seq_files(find(contains(seq_names, seq_name)));
        seq_name = file.name;
        load(fullfile(file.folder, seq_name))
        seq_name = extractBefore(seq_name,'.mat');

        % make new signals structure
        new_sig_struct = signals.(seq_name);

        % collect indices where to insert b0
        pos = [];

        for i = unique(seq.gc_ind)
            for j = unique(seq.tm_ind)
                pos = [pos min(find(seq.gc_ind == i & seq.tm_ind == j))];
            end
        end



        tmp_seq.gc_ind   = insert_copy2array(seq.gc_ind, pos);
        tmp_seq.tm_ind   = insert_copy2array(seq.tm_ind, pos);
        tmp_seq.tm       = insert_copy2array(seq.tm, pos);


        tmp_seq.gc       = insert_copy2array(seq.gc, pos);
        tmp_seq.qc       = insert_copy2array(seq.qc, pos);


        tmp_seq.isDDE    = insert_copy2array(seq.isDDE, pos);

        tmp_seq.b        = insert_2array(seq.b, pos*0, pos);
        tmp_seq.b1       = insert_2array(seq.b1, pos*0, pos);
        tmp_seq.b2       = insert_2array(seq.b2, pos*0, pos);

        tmp_seq.bm       = insert_2array(seq.bm, pos*0, pos);
        tmp_seq.g1       = insert_2array(seq.g1, pos*0, pos);
        tmp_seq.g2       = insert_2array(seq.g2, pos*0, pos);

        tmp_seq.delta    = insert_copy2array(seq.delta, pos);

        tmp_seq.n = length(tmp_seq.b);

        tmp_seq.b_array = unique(tmp_seq.b);
        tmp_seq.b_ind = zeros(1,tmp_seq.n);
        tmp_seq.nb = length(tmp_seq.b_array);
        for i = 1:tmp_seq.nb
            tmp_seq.b_ind(tmp_seq.b == tmp_seq.b_array(i)) = i;
        end


        % copy missing fields
        f = fields(seq);
        for i = 1:numel(f)
            if ~isfield(tmp_seq,f{i})
                tmp_seq.(f{i}) = seq.(f{i});
            end
        end


        seq = tmp_seq;
        % save new sequence
        new_seq_name = [seq_name '_b0'];
        save(fullfile(file.folder, new_seq_name),'seq')

        % get all relevant fields
        f = fields(new_sig_struct);
        f = f(contains(f,'_k'));
        N = numel(f);

        for n = 1:N;
            s = signals.(seq_name).(f{n});
            new_signals.(seq_name).(f{n}) = insert_2array(s, ones(length(pos),1), pos);
        end

    end

    new_signals = renameStructField(new_signals, seq_name, new_seq_name);

    % save new signals
    signals = new_signals;
    new_signal_fn = strrep(signal_fn,'.mat','_b0.mat');
    save(new_signal_fn,'signals');

end

function new_array = insert_2array(array, val, pos)
new_array = [];

for i = 1:length(array)
    pos_ind = find(pos == i);

    if isempty(pos_ind)
        new_array = [new_array array(i)];
    else
        new_array = [new_array val(pos_ind) array(i)];
    end
end
end

function new_array = insert_copy2array(array, pos)
new_array = [];

for i = 1:length(array)
    pos_ind = find(pos == i);

    if isempty(pos_ind)
        new_array = [new_array array(i)];
    else
        new_array = [new_array array(i) array(i)];
    end
end
end



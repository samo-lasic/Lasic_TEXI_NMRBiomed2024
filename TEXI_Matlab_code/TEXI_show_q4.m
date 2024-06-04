% help script to check q4 waveforms

clear all
close all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

seq_folder = fullfile(root_dir, 'seq_TEXI');

fn_append = '';  % ignore if empty
select_seq_names = {'FEXI_preclin7_dz500_tm400_bf2000_b1300','VT_preclin7_dz1000_tm250_b2600'};
select_seq_names = select_seq_names(1);

filter_names = {}; % ignore if empty

tm_ind = [2]; % if empty do all
gc_ind = []; % if empty do all
b_ind = []; % if empty do max b

dt = 1e-6;
N0pad = 1e7;

opt.fig_n = 3; % ofset figure number
opt.fig_width_mm = 40;
opt.fig_aspect = .4;
opt.display_ratio = .15;

opt.fLIM = [0 100];
opt.lw = 3;
opt.ms = 70;
opt.fs = 24;
opt.title_fs_ratio = .7;
opt.show_all = 0;
opt.show_labels = 1;

files = dir(seq_folder);
ind = ones(1, numel(files));

if ~isempty(fn_append)
    ind = ind & contains({files.name}, fn_append); 
end

% optionally add seq name filter here
if ~isempty(select_seq_names)
    ind = ind & contains({files.name}, select_seq_names);

    if ~isempty(filter_names)
        ind = ind & ~contains({files.name}, filter_names);
    end
end

ind = find(ind & ~contains({files.name},{'_DvsR','b_'}));
files = files(ind);

for nSeq = 1:numel(files)
    seq_file = files(nSeq);
    seq_path = fullfile(seq_file.folder, seq_file.name);


    wfm = make_wfm(seq_path, dt);

    load(seq_path)

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

    sub_ind = find(sub_ind);
    sub_seq = get_sub_seq_from_ind(seq, sub_ind);

    %     Tenc_max = 2*(max(sub_seq.delta) + sub_seq.delta_pi + sub_seq.delta_max + sub_seq.delta_c) + max(sub_seq.tm);

    opt.size_mm = opt.fig_width_mm * [1 opt.fig_aspect * sub_seq.n];

    % find max length
    L = 0;
    for c = 1:sub_seq.n
        l = length(wfm.seq(sub_ind(c)).q);
        if l > L
            L = l;
        end
    end

    q4_array = zeros(L,sub_seq.n);

    figure(1), clf, hold on
    % fill up q4_array
    for c = 1:sub_seq.n
        q = wfm.seq(sub_ind(c)).q;
        
        plot(q)

        q2 = q.^2;
        l = length(q2);
        q4 = xcorr(q2,q2) * dt;
        q4 = [0; q4((l+1):end)];

        % normalize
        b = sum(q2) * dt;
        q4_array(1:l,c) = q4/b^2;
        %         sum(q4/b^2)*dt
    end

    ub_ind = unique(sub_seq.b_ind);
    Nb = length(ub_ind);
    n_tm_gc = sub_seq.n/Nb;

    % b-color
    bc = linspace(0.2,1,Nb);


    % normalize
    % q4_array = q4_array./max(q4_array(:));


    % plot q4 waveforms

    % color_order = [1 0 0; 0 0 1; .5 0 0 ; 0 0 .5];
    color_order = turbo(n_tm_gc);

    t = 1e3*linspace(0,L*dt,L);
    fh = figure(length(findobj('type','figure')) + opt.fig_n + 1);
    %fh = figure;
    clf
    set(gcf, 'InvertHardCopy', 'off')
    set(fh,'Color','white', 'Units', 'inches', ...
        'PaperPosition',[0 0 opt.size_mm / 25.4 ], 'PaperPositionMode', 'auto');
    screen_size = get(0,'screensize');
    fh.Position(3:4) = opt.display_ratio * opt.size_mm;

    for c = 1:n_tm_gc

        subplot(n_tm_gc,1,c)
        hold on
        for nb = 1:Nb
            tm_gc_ind = find(sub_seq.b_ind == ub_ind(nb));

            if nb == 1
                title_str = sprintf('tm = %d ms qc = %g m^{-1}', round(sub_seq.tm(tm_gc_ind(c)) * 1e3), sub_seq.qc(tm_gc_ind(c)));
            end

            plot(t,q4_array(:,tm_gc_ind(c)),'-k','LineWidth',opt.lw,'MarkerSize',opt.ms,'color',bc(nb) * color_order(c,:));
        end

        xlim(1e3*[0*dt dt*L])

        th = title(title_str);

        if c == 1
            ylabel('q4(t)/b^2 [s^-1]')
        elseif c == sub_seq.n
            xlabel('time [ms]')
        end
        set(gca,'YScale', 'lin','LineWidth',opt.lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',opt.fs)
        th.FontSize = th.FontSize * opt.title_fs_ratio;

    end

end


clear all
close all

% Make parameters for FEXI and TEXI sequences. 
% 1 - make sequence with b*Vw = const. and calculate/save ADCs
% 2 - make sequence with b = const. (match max delta) and calculate/save ADCs (not used in the manuscript, but could be used for different tunings not based on low-freq. approx.)
% 3 - tune sequence with b = const. based on the ADCs from the previous step (not used in the manuscript, but could be used for different tunings not based on low-freq. approx.)
% 4 - calculate/save ADCs for the tuned sequence

% output:
% VT_*.mat
% b_*.mat
% ST_*.mat
% FEXI_*.mat

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));

gmr = 26.75e7;

fn_prepend = '';
seq_folder = fullfile(root_dir, 'seq_TEXI');

% ------------------
steps = 4 %[1:4]; % which steps to run (steps 1 and 4 are used in the paper)
do_save_ADCs = 0;

fn_append = 'preclin7';
% fn_append = 'preclinX';
constraints.g_lim = 600e-3;
constraints.delta_min = 3e-3;
constraints.delta_pi = 3e-3;

% ------------------


% waveform time resolustion (only used for ADC calculations and checking waveforms)
dt = 1e-6;

% --------------- parameters used for tuning -------------
D0 = 2e-9;
R_array = [.5 2.5 5] * 1e-6; % this is used to generate signals based on the GPA models (not MC), can be set to values used for tuning, but can be different from R_tune_array
R_tune_array = R_array(2);
restriction_geometry = 1; %1-spherical, 2-cylindrical, 3-planar

N0pad = 1e7; % zero padding FFT

% this was used in some tests, setting it to 1 may not work
invert_polarity = 0; % 0 or 1 to have positive or negative polarity of the detection block


for step = steps
    switch step
        case 1

            g_lim = constraints.g_lim; delta_pi = constraints.delta_pi;
            delta_lim = [constraints.delta_min 50e-3]; % upper limit is just for showing extended range in finding suitable Vw

            dz_min = 1 * 1e-3;
            tm_array = round(linspace(5,250,5))*1e-3;
            b_array = linspace(300, 2600, 7)*1e6;


            % --------------------------------------------------
            inlude_DDE = 1; % use one DDE (only gradients in the first block + mixing), other are TDE


            % ------- crusher calculation (with slice) ----------

            % delta_c = max(qc_array)/gmr/g_lim; % using minimum possible delta_c
            delta_c = 2e-3;


            dfds = 6; % ﻿Lasič et al., Magn Reson Med. 2018;79(4):2228–35.

            dz = [Inf 5 2.5 1] * 1e-3;
            [qc_array, gc_array] = slice2butterfly(dz, delta_c, dfds);

            %         gc_array = qc_array/gmr/delta_c;
            % --------------------------------------

            %                 we get the same qb if delta_s = 2e-3, df = 3000 and delta_c = 1e-3 or if delta_c_effective = 1.6e-3 and heard 90 pulse is used
            %                 qb = gmr*(max(gc_array)*delta_c + gs*delta_s/2)
            %                 max(gc_array)*delta_c + gs*delta_s/2 = max(gc_array)*delta_c_effective
            %                 delta_c_effective = (max(gc_array)*delta_c + gs*delta_s/2) / max(gc_array);
            %                 delta_c_effective = 1.6e-3
            %                 gmr*max(gc_array)*delta_c_effective
            %                 max(qc_array)



            % ------------  find suitable Vw --------------

            [Vw, b_array] = find_suitable_Vw_and_b(gc_array, qc_array, tm_array, b_array, delta_c, delta_pi, inlude_DDE, g_lim, delta_lim);


            % ----------------------------  sequences with constant Vw -------------------

            % ---------------------- save sequences ---------------------------
            if ~isnan(Vw)
                seq_name = sprintf('%s_dz%i_tm%d_b%d', fn_append, round(dz_min*1e6), round(max(tm_array)*1e3), round(max(b_array)*1e-6));
                seq_name = make_seq_name(fn_prepend, 'VT', seq_name);

                do_save_seq = 1;
                seq = save_seq_with_const_Vw(seq_folder, seq_name, gc_array, delta_c, delta_pi, b_array, tm_array, Vw, inlude_DDE, do_save_seq);
                % ---------------------- save ADCs ---------------------------
                if do_save_ADCs
                    save_ADCs(seq_folder, seq_name, invert_polarity, D0, R_array, N0pad, dt, restriction_geometry)
                end

            end


            % ------------------------ save tuned sequences based on V(omega) --------------------------------------------------
            % use delta_max and other parameters such as delta_c and delta_pi from the saved VT sequences

        case 2
            % ---------------------- save sequences with constant b along with the ADCs ---------------------------
            % warning: calculating ADC table may take a while

            seqName_in = 'VT';
            seqName_out = 'b';
            Ndelta = 30; % delta grid used for tuning
            inlude_DDE = 1;
            show_b = 0;


            files = dir(seq_folder);
            ind = contains({files.name},sprintf('%s_%s_',seqName_in, fn_append));

            % optionally add seq name filter here
            %     ind = ind & contains({files.name},{'VT_preclin5_dz500_tm400_b2600','VT_preclin5_dz500_tm200_b2600'});

            ind = find(ind & ~contains({files.name},{'_DvsR','ST_','b_'}));
            files = files(ind);

            for c = 1:numel(files)
                fn_ref_seq = fullfile(files(c).folder,files(c).name);

                [path, fn_out_seq] = fileparts(fn_ref_seq);
                fn_out_seq = ['b_' extractAfter(fn_out_seq,'_')];
                fn_out_seq = fullfile(path, fn_out_seq);

                save_seq_with_const_b(fn_ref_seq, fn_out_seq, constraints, inlude_DDE, Ndelta, dt, seqName_out, show_b);

                % ---------------------- save ADCs for sequences with constant b ---------------------------
                % ADCs for all the sizes in R_tune_array
                if do_save_ADCs
                    [~, seq_name] = fileparts(fn_out_seq);
                    save_ADCs(seq_folder, seq_name, invert_polarity, D0, R_tune_array, N0pad, dt, restriction_geometry)
                end


            end


        case 3
            % ---------------------- size tune sequences ---------------------------

            seqName_in = 'b';


            files = dir(seq_folder);
            ind = contains({files.name}, sprintf('%s_%s_', seqName_in, fn_append));

            % optionally add seq name filter here
            %ind = ind & contains({files.name},{'b_preclin5_dz500_tm400_b2600','b_preclin5_dz500_tm200_b2600'});

            ind = find(ind & ~contains({files.name},{'_DvsR','VT_','ST_'}));
            files = files(ind);

            tol = 1; % relative ADC difference (percent) from target in find_tunable
            tune_step = 1e-3; % for iterating refine_tuneADC_array

            % check tuning process with plots
            show_ADC = 0;
            show_refine_steps = 0;
            show_tunable_in_full_range = 0; % just check full range n_tunable

            for c = 1:numel(files)
                fn_seq2tune = fullfile(files(c).folder,files(c).name);

                % just checking
                % load(fn_seq2tune), display_seq_par_range(files(c).name, seq, 0);


                for n = 1:length(R_tune_array)
                    R_tune = R_tune_array(n);
                    fn_out_seq = size_tune_seq(R_tune, fn_seq2tune, tol, constraints, tune_step, show_ADC, show_refine_steps, show_tunable_in_full_range);


                    % ---------------------- save ADCs ---------------------------
                    if do_save_ADCs
                        [~, seq_name] = fileparts(fn_out_seq);

                        % just checking
                        % load(fn_out_seq), display_seq_par_range(seq_name, seq, 0);

                        save_ADCs(seq_folder, seq_name, invert_polarity, D0, R_array, N0pad, dt, restriction_geometry)
                    end

                end

            end



        case 4
            % ------------------------ FEXI -----------------------
            % warning: check manually that max gradient is not exceeded


            fn_append = 'preclin7'; % the rest of string is appended

            g_lim = constraints.g_lim;
            delta_pi = constraints.delta_pi;
            dz_min = .5 * 1e-3;

            tm_array = round(linspace(5,400,5))*1e-3;
            delta = 5e-3;
            bf = 2000*1e6;
            b_array = linspace(50e6, 1300*1e6, 4);


            % -------- general ---------------------
            Delta = delta_pi + delta;
            delta_max = delta;

            % ------- crusher calculation (with slice) ----------

            % delta_c = max(qc_array)/gmr/g_lim;

            delta_c = 2e-3;
            Nc = 5;
            dz_max = 5 * 1e-3; % only used when Nc > 2
            dfds = 6; % ﻿Lasič et al., Magn Reson Med. 2018;79(4):2228–35.
            [qc_array, gc_array, dz] = make_crusher_array(Nc, dfds, delta_c, dz_min, dz_max);
            qc_array = fliplr(qc_array);
            gc_array = fliplr(gc_array);


            if (1)
                fn_append = 'preclinX'; % the rest of string is appended

                g_lim = constraints.g_lim;
                delta_pi = constraints.delta_pi;
                dz_min = 1 * 1e-3;

                tm_array = [5 100 250 400]*1e-3;
                delta = 5e-3;
                bf = 2000*1e6;
                b_array = linspace(50e6, 1300*1e6, 4);


                % -------- general ---------------------
                Delta = delta_pi + delta;
                delta_max = delta;

                % ------- crusher calculation (with slice) ----------

                % delta_c = max(qc_array)/gmr/g_lim;

                delta_c = 2e-3;

                dfds = 6; % ﻿Lasič et al., Magn Reson Med. 2018;79(4):2228–35.

                dz = [Inf 5 2.5 1] * 1e-3;
                [qc_array, gc_array] = slice2butterfly(dz, delta_c, dfds);
            end


            % --------------------------------------

            fn_append = sprintf('%s_dz%i_tm%d_bf%d_b%d', fn_append, round(dz_min*1e6), round(max(tm_array)*1e3), round(bf*1e-6), round(max(b_array)*1e-6));


            seq = make_seq_FEXI_gc(gc_array, delta_c, delta_pi, delta_max, delta, bf, b_array, tm_array);

            if (0) % checking sequence
                dt = 1e-6;
                for c = 1:seq.n

                    [g, q] = def_FEXI(seq.delta, seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);

                    figure(1),clf
                    subplot(2,1,1), plot(g)
                    subplot(2,1,2), plot(q)

                    q2 = q.^2;

                    b(c) = sum(q2)*dt;
                end

                figure(1),clf
                hold on
                plot(seq.b,'o')
                plot(b,'.')
            end

            fn_out_seq = make_seq_name(fn_prepend, 'FEXI', fn_append);

            display_seq_par_range(fn_out_seq, seq, 1)

            mkdir(seq_folder)
            save(fullfile(seq_folder, fn_out_seq),'seq')

            % ---------------------- save ADCs ---------------------------

            if do_save_ADCs
                [~, seq_name] = fileparts(fn_out_seq);
                save_ADCs(seq_folder, seq_name, invert_polarity, D0, R_array, N0pad, dt, restriction_geometry)
            end
    end
end
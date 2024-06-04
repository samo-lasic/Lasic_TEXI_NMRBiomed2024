% Simulate blood-brain barrier exchange with Gaussian IVIM diffusion contrast in FEXI with crushers. 
% Exchange rates are estimated with the original AXR fit and with the Ning model.

clear all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));


simulate = 0; % if 0 just load

sim_cases = [9];


for sim_case = sim_cases
    switch sim_case
        case 1
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 300 * 1e-3, 300 * 1e6, 50 * 1e6, 1300 * 1e6, 1.5 * 1e-3, 5 * 1e-3);

        case 2
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 500 * 1e-3, 300 * 1e6, 50 * 1e6, 1300 * 1e6, 1.5 * 1e-3, 5 * 1e-3);

        case 3
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 300 * 1e-3, 300 * 1e6, 50 * 1e6, 2300 * 1e6, 1.5 * 1e-3, 5 * 1e-3);

        case 4
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 500 * 1e-3, 300 * 1e6, 50 * 1e6, 2300 * 1e6, 1.5 * 1e-3, 5 * 1e-3);


        case 5
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 300 * 1e-3, 300 * 1e6, 50 * 1e6, 1300 * 1e6, 1 * 1e-3, 5 * 1e-3);

        case 6
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 500 * 1e-3, 300 * 1e6, 50 * 1e6, 1300 * 1e6, 1 * 1e-3, 5 * 1e-3);

        case 7
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 300 * 1e-3, 300 * 1e6, 50 * 1e6, 2300 * 1e6, 1 * 1e-3, 5 * 1e-3);

        case 8
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 500 * 1e-3, 300 * 1e6, 50 * 1e6, 2300 * 1e6, 1 * 1e-3, 5 * 1e-3);


        case 9
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 300 * 1e-3, 250 * 1e6, 50 * 1e6, 500 * 1e6, 1 * 1e-3, 5 * 1e-3);

        case 10
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 500 * 1e-3, 250 * 1e6, 50 * 1e6, 500 * 1e6, 1 * 1e-3, 5 * 1e-3);

        case 11
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 300 * 1e-3, 250 * 1e6, 50 * 1e6, 1000 * 1e6, 1 * 1e-3, 5 * 1e-3);

        case 12
            [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
                set_par(6.5 * 1e-9, 0.65 * 1e-9, 0.05, 500 * 1e-3, 250 * 1e6, 50 * 1e6, 1000 * 1e6, 1 * 1e-3, 5 * 1e-3);


    end


    fn = sprintf('sim_%d_FEXI_IVIM',sim_case);


    if simulate

        k_array = [0:2:20]; %[0:2:20];

        D = [D_intra 0 ; 0 D_extra];
        f = [f_intra 1-f_intra]';


        g_lim = 0.6;

        delta_pi = 3e-3;

        tm_array = round(linspace(5,tm_max*1e3,5))*1e-3;
        delta = 5e-3;
        b_array = linspace(b_min, b_max, 4); %linspace(50e6, 1300*1e6, 4);


        Delta = delta_pi + delta;
        delta_max = delta;

        % ------- crusher calculation (with slice) ----------


        delta_c = 2e-3;
        Nc = 5;
        dfds = 6; % ﻿Lasič et al., Magn Reson Med. 2018;79(4):2228–35.
        [qc_array, gc_array, dz] = make_crusher_array(Nc, dfds, delta_c, dz_min, dz_max);
        qc_array = fliplr(qc_array);
        gc_array = fliplr(gc_array);

        seq = make_seq_FEXI_gc(gc_array, delta_c, delta_pi, delta_max, delta, b_f, b_array, tm_array);
        

        dt = 1e-5;


        Nk = length(k_array);

        sig_array = zeros(Nk,seq.n);
        tic
        parfor c = 1:seq.n

            [g, q] = def_FEXI(seq.delta, seq.delta_max, seq.g1(c), seq.g2(c), seq.tm(c), seq.gc(c), seq.delta_c, seq.delta_pi, dt);

            q2 = q.^2;


            for nk = 1:Nk

                k = k_array(nk);
                kie = k * (1-f_intra);
                kei = k * f_intra;
                K = [-kie kei ; kie -kei];

                S = f;
                for nt = 1:length(q)
                    S = expm((-q2(nt)*D + K)*dt)*S;
                end
                sig_array(nk,c) = sum(S);
            end
        end
        toc

        sim.D_intra = D_intra;
        sim.D_extra = D_extra;
        sim.f_intra = f_intra;
        sim.k_array = k_array;
        sim.dz = dz;
        sim.gc_array = gc_array;
        sim.qc_array = qc_array;
        sim.Nc = length(qc_array);
        sim.Nk = Nk;
        sim.dt = dt;
        sim.seq = seq;
        sim.sig_array = sig_array;


        save(fn,'sim')

    else

        load(fn,'sim')
    end % --------------------------------------------------------------


    fit_model_array = [1 2]; %[1 2]

    do_plot_loop = 0;
    do_save_fit_sig = 1;
    n_inits = 10;


    seq = sim.seq;
    sig_array = sim.sig_array;
    Nc = sim.Nc;
    Nk = sim.Nk;


    for fit_model = fit_model_array


        if fit_model == 1
            fitP_array = zeros(Nc, Nk, 1, 4 + length(seq.tm_array));
        elseif fit_model == 2
            fitP_array = zeros(Nc, Nk, 1, 3);
        elseif fit_model == 3
            fitP_array = zeros(Nc, Nk, 1, 4);
        end

        parfor nc = 1:Nc


            [sub_ind, sub_seq] = get_sub_seq_gc_b2(seq, nc, 0, Inf);

            sub_sig = sig_array(:, sub_ind);

            if fit_model == 1 % FEXI
                xps = make_FEXI_xps(sub_seq);
                opt = fexi11_opt();
                opt.fexi11.do_plot = do_plot_loop;

            elseif fit_model == 2 %  Ning_h_norm (fix S0)
                opt = Ning_h_norm_opt();
                opt.Ning_h_norm.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, 1, sim.dt); % exchange sensitivity

            elseif fit_model == 3
                xps = make_FEXI_xps(sub_seq);
                aux = expm_aux(sub_seq, 1, sim.dt); % exchange sensitivity
                opt = expm_opt();
                opt.expm.do_plot = do_plot_loop;

            end

            opt.n_inits = n_inits; % number of initial guesses

            for nk = 1:Nk
                sig = squeeze(sub_sig(nk, :))';

                if fit_model == 1 % FEXI
                    p = fexi11_data2fit(sig, xps, opt);

                    if do_save_fit_sig || do_plot_loop
                        sFIT = fexi11_fit2data(p, xps);
                    end
                    if do_plot_loop
                        fixi_fig(1, 1, sub_seq, sig, sFIT, p);
                        pause(.1)
                    end

                elseif fit_model == 2  % Ning_h_norm (fix S0)
                    p = Ning_h_norm_data2fit(sig', aux, opt);

                    if do_save_fit_sig || do_plot_loop
                        sFIT = Ning_h_norm_fit2data(p,aux);
                    end
                    if do_plot_loop
                        fixi_fig(2, 1, sub_seq, sig, sFIT, p);
                        pause(.1)
                    end
                elseif fit_model == 3  % expm (fix S0)
                    tic
                    sFIT1 = expm_fit2data([1 0.9 1e-9 10e-9],aux);
                    toc
                    p = expm_data2fit(sig', aux, opt);
                    if do_save_fit_sig || do_plot_loop
                        sFIT = expm_fit2data(p,aux);
                    end
                    if do_plot_loop
                        fixi_fig(2, 1, sub_seq, sig, sFIT, p);
                        pause(.1)
                    end


                end


                fitP_array(nc,nk,1,:) = p;

                fit_sig_struct(nc,nk,1).sFIT = sFIT(:);
                fit_sig_struct(nc,nk,1).sub_ind = sub_ind;
                % display(sprintf('nc/Nc = %.2f, nk/Nk = %.2f, nR/NR = %.2f',nc/Nc, nk/Nk,nR/NR))

            end
        end


        % optionally save fit signal
        if do_save_fit_sig
            fit_sig_array = 0*sig_array;
            for nc = 1:Nc
                for nk = 1:Nk
                    sub_ind = fit_sig_struct(nc,nk,1).sub_ind;
                    sFIT = fit_sig_struct(nc,nk,1).sFIT;
                    fit_sig_array(nk, sub_ind) = sFIT;

                    sig = squeeze(sig_array(nk,sub_ind))';
                    sum_square_res(nc,nk,1) = sum((sig-sFIT).^2);
                end
            end

            res.sum_square_res = sum_square_res;
            res.sig_array = reshape(sig_array,size(sig_array,1),1,size(sig_array,2));
            res.fit_sig_array = reshape(fit_sig_array,size(fit_sig_array,1),1,size(fit_sig_array,2));
        end

        res.in.Nc = Nc;
        res.in.Nk = Nk;
        res.in.gc_array = seq.gc_array;
        res.in.qc_array = seq.qc_array;
        res.in.k_array = sim.k_array;

        res.in.dt = sim.dt;

        res.fitP_array = fitP_array;

        res_fn = strcat(fn, '_fit',num2str(fit_model));

        res.in.fn = res_fn;

        path = fullfile('res_IVIM');
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        mkdir(path)
        res_path = fullfile(path, res_fn);
        save(res_path,'res')
        display(sprintf('%s',res_path))

    end
end

function [D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max] = ...
    set_par(D_intra, D_extra, f_intra, tm_max, b_f, b_min, b_max, dz_min, dz_max)
end


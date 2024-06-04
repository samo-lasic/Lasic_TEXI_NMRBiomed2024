% Estimate fit precision for selected data. 

% note on subsampling simulation data:
% this is a workaround, which could be adressed by setting more specific names (with bmin info)
% 1. make simulation data and VT sequences with increased lower b
% use: TEXI_make_data_and_seq_reduce_b_VT.m
% save these data in a separate folder
% save also sub-sequence for later use
% "reduceB" is appended to file names
% 2. fit data
% make sure to use the correct data with "reduceB" appended
% make sure to use the correct sequence with "reduceB" appended
% 3. results
% save results with "reduceB" appended to a separate folder
% 4. plot
% make sure to use the correct results with "reduceB" appened



clear all
close all


warning('off','all')
restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

root_dir = fileparts(fileparts(pwd));


fn_prepend = '';
fn_append = 'clin7';
fn_substrate_tag = '07';

simulation_data_folder = fullfile('data_TEXI','simulation');
seq_folder = 'seq_TEXI';
results_folder = fullfile('data_TEXI','results',['sim' fn_substrate_tag]);
%results_folder = fullfile('results',['sim' fn_substrate_tag]);

% % VT only
% sim_fn_array = {...
%     'signals_spheres_d1_regular_rho50',...
%     'signals_spheres_d3_regular_rho50',...
%     'signals_spheres_d5_regular_rho50',...
%     'signals_spheres_d7_regular_rho50',...
%     'signals_spheres_d9_regular_rho50',...
%     'signals_spheres_d11_regular_rho50'...
%     };

% data_name = 'spheres_d3_VT';
% 
% protocol = {};
% protocol(end+1).sim_fn = 'signals_spheres_d3_regular_rho50';
% protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600';
% protocol(end).b_max = Inf;
% 
% protocol(end+1).sim_fn = 'signals_spheres_d3_regular_rho50';
% protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600';
% protocol(end).b_max = 1450;
% 
% protocol(end+1).sim_fn = 'signals_spheres_d3_regular_rho50_bmin1450';
% protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600_bmin1450';
% protocol(end).b_max = 1833;
% 
% protocol(end+1).sim_fn = 'signals_spheres_d3_regular_rho50_bmin1834';
% protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600_bmin1834';
% protocol(end).b_max = 1450;


d_str = 'd11'

data_name = sprintf('spheres_%s_VT',d_str);

protocol = {};
protocol(end+1).sim_fn = sprintf('signals_spheres_%s_regular_rho50',d_str);
protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600';
protocol(end).b_max = Inf;

protocol(end+1).sim_fn = sprintf('signals_spheres_%s_regular_rho50',d_str);
protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600';
protocol(end).b_max = 1450;

protocol(end+1).sim_fn = sprintf('signals_spheres_%s_regular_rho50_bmin1450',d_str);
protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600_bmin1450';
protocol(end).b_max = 1833;

protocol(end+1).sim_fn = sprintf('signals_spheres_%s_regular_rho50_bmin1834',d_str);
protocol(end).seq_name = 'VT_preclin7_dz1000_tm250_b2600_bmin1834';
protocol(end).b_max = 1450;


save_data = 1; % if 0 load

k_sel = 10;

% use Inf for no limit, only applies to Ning fit of non-FEXI data
% b_max_array = 50 + [1450   1833    2217    2600];
%b_max_array = 50 + [1450   1833];
%b_limit_array = Inf;

% b_min = 0*300; % test using only  b-values > b_min for FEXI detection (b2)

% 1-FEXI, 2-Ning_h, 3-Ning_h_norm, 4-Ning_Gamma, 5-Ning_g,
% 6-Ning_h3_norm, 7-Ning_g_norm, 8-Ning_h3b_norm, 9-Ning_h_intra_norm

fit_model_for_FEXI_data = 3; %[1 3]; %[1 3 7];
fit_model_for_nonFEXI_data = 3; %[3 7];

dt = 1e-5; % waveform resolution used in simulations

do_plot_loop = 0; % for debugging

do_save_fit_sig = 0;

n_inits = 10; %10; % number of inital guesses
SNR = 100;
N_SNR = 1000; %1000

data_path = sprintf('fitP_SNR%d_N%d_inits%d_%s',SNR,N_SNR,n_inits,data_name);
data_path = fullfile(fileparts(fileparts(pwd)), results_folder, data_path);
display(sprintf('%s',data_path))

if save_data
    fitP_data = {};

    Np = numel(protocol);
    for n_protocol = 1:Np
        tic

        sim_fn = protocol(n_protocol).sim_fn;
        seq_name = protocol(n_protocol).seq_name;
        b_max = protocol(n_protocol).b_max;

        % load simulated signals
        signal_fn = fullfile(root_dir, simulation_data_folder, fn_substrate_tag, char(sim_fn));
        load(signal_fn)

        % load sequence
        seq_fn = fullfile(root_dir, seq_folder,seq_name);
        load(seq_fn)

        isFEXI = contains(seq_name,{'FEXI'});


        % make arrays for fitting
        if isFEXI % b < bmax for non-FEXI sequences
            b_max = Inf;
            fit_model = fit_model_for_FEXI_data;
        else
            fit_model = fit_model_for_nonFEXI_data;
        end


        % make arrays from simulation output structure
        [sig_array, d_array, k_array, NR, Nk] = sim2array(signals, seq_name);

        % select k
        if ~isempty(k_sel)
            [~,k_ind] = min(abs(k_array-k_sel));
            k_array = k_array(k_ind);
            Nk = 1;
            sig_array = sig_array(k_ind,:,:);
        end

        R_array = 0.5*1e-6 * d_array;

        Nc = seq.nc;

        % ----------  fitting  ----------------------------

        % 1-FEXI, 2-Ning_h, 3-Ning_h_norm, 4-Ning_Gamma, 5-Ning_g,
        % 6-Ning_h3_norm, 7-Ning_g_norm, 8-Ning_h3b_norm, 9-Ning_h_intra_norm

        if fit_model == 1
            fitP_array = zeros(Nc, Nk, NR, 4 + length(seq.tm_array)); % + non-filtered
        elseif fit_model == 2 || fit_model == 4 || fit_model == 5
            fitP_array = zeros(Nc, Nk, NR, 3 + length(seq.tm_array));
        elseif fit_model == 3 || fit_model == 7
            fitP_array = zeros(Nc, Nk, NR, 3);
            fitP_std_array = zeros(Nc, Nk, NR, 3);

        elseif fit_model == 6 || fit_model == 8 || fit_model == 9
            fitP_array = zeros(Nc, Nk, NR, 4);
        else
            display('ERROR ... check fitting model')
            return
        end

        parfor nc = 1:Nc
%         for nc = 1:Nc
            % display(sprintf('nc/Nc = %.2f',nc/Nc))
            if isFEXI
                [sub_ind, sub_seq] = get_sub_seq_gc_b2(seq, nc, b_min, Inf);
            else
                [sub_ind, sub_seq] = get_sub_seq_gc_b(seq, nc, 0, b_max);
            end

            sub_sig = sig_array(:, :, sub_ind);


            if fit_model == 1 % FEXI
                xps = make_FEXI_xps(sub_seq);
                opt = fexi11_opt();
                opt.fexi11.do_plot = do_plot_loop;
            elseif fit_model == 2 % Ning_h
                opt = Ning_h_opt();
                opt.Ning_h.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 3 %  Ning_h_norm (fix S0)
                opt = Ning_h_norm_opt();
                opt.Ning_h_norm.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 4 % Ning_Gamma (approximation for low exchange weighting)
                opt = Ning_Gamma_opt();
                opt.Ning_Gamma.do_plot = do_plot_loop;
                aux = Ning_Gamma_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 5 %  Ning_g (with gamma distribution)
                opt = Ning_g_opt();
                opt.Ning_g.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 6 % Ning_h3_norm (S0 = 1)
                opt = Ning_h3_norm_opt();
                opt.Ning_h3_norm.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 7 %  Ning_g_norm (with Gamma distribution, S0 = 1)
                opt = Ning_g_norm_opt();
                opt.Ning_g_norm.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 8 % Ning_h3b_norm (with 3 cumulants + scalling h with b, S0 = 1)
                opt = Ning_h3b_norm_opt();
                opt.Ning_h3b_norm.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            elseif fit_model == 9 % Ning_h3_intra_norm (fix S0 + offset due to the intra-compartmental variance of substrate mean diffusivities, used for directional/powder averaging)
                opt = Ning_h_intra_norm_opt();
                opt.Ning_h_intra_norm_opt.do_plot = do_plot_loop;
                aux = Ning_h_aux(sub_seq, isFEXI, dt); % exchange sensitivity
            end



            opt.n_inits = n_inits; % number of initial guesses


            for nk = 1:Nk
                for nR = 1:NR
                    sig = squeeze(sub_sig(nk, nR, :));

                    if fit_model == 1 % FEXI
                        p = fexi11_data2fit(sig, xps, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = fexi11_fit2data(p, xps);
                        end
                        if do_plot_loop
                            fixi_fig(1, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end


                    elseif fit_model == 2 % Ning_h

                        p = Ning_h_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_h_fit2data(p,aux);
                        end

                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end

                    elseif fit_model == 3 % Ning_h_norm (fix S0)
                        %                     p = Ning_h_norm_data2fit(sig', aux, opt);

                        % SNR
                        sIN = sig';
                        p_SNR = zeros(N_SNR,3);
                        for nSNR = 1:N_SNR
                            s2fit = abs( sIN + 1/SNR * randn(size(sIN)) );
                            p_SNR(nSNR,:) = Ning_h_norm_data2fit(s2fit, aux, opt);
                        end
                        p = mean(p_SNR);
                        p_std = std(p_SNR);


                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_h_norm_fit2data(p,aux);
                        end
                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end

                    elseif fit_model == 4 % Ning_Gamma
                        p = Ning_Gamma_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_Gamma_fit2data(p,aux);
                        end
                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end

                    elseif fit_model == 5 % Ning_g (gamma distribution)

                        p = Ning_g_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_g_fit2data(p,aux);
                        end

                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end



                    elseif fit_model == 6 % Ning_h3_norm

                        p = Ning_h3_norm_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_h3_norm_fit2data(p,aux);
                        end

                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            display(sprintf('c2 = %g, c3 = %g, c3/c2 = %g', p(2), p(3), p(3)/p(2)))
                            pause(.1)
                        end


                    elseif fit_model == 7 % Ning_g_norm (gamma distribution)

                        p = Ning_g_norm_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_g_norm_fit2data(p,aux);
                        end

                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end

                    elseif fit_model == 8 % Ning_h3b_norm (with 3 cumulants + scalling h with b, S0 = 1)

                        p = Ning_h3b_norm_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_h3b_norm_fit2data(p,aux);
                        end

                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            display(sprintf('c2 = %g, c3 = %g, c3/c2 = %g', p(2), p(3), p(3)/p(2)))
                            pause(.1)
                        end

                    elseif fit_model == 9 % Ning_h_intra_norm (fix S0)
                        p = Ning_h_intra_norm_data2fit(sig', aux, opt);

                        if do_save_fit_sig || do_plot_loop
                            sFIT = Ning_h_intra_norm_fit2data(p,aux);
                        end
                        if do_plot_loop
                            fixi_fig(2, isFEXI, sub_seq, sig, sFIT, p);
                            pause(.1)
                        end

                    end

                    if do_plot_loop
                        display(sprintf('nc/Nc = %d/%d, nk/Nk = %d/%d, nR/NR = %d/%d', ...
                            nc, Nc, nk, Nk, nR, NR))
                        %pause(1);
                    end

                    fitP_array(nc,nk,nR,:) = p;
                    fitP_std_array(nc,nk,nR,:) = p_std;


                end


            end

        end


        time = toc;
        display(sprintf('protocol %d/%d done in %.1f min (%.1f h)', n_protocol, Np, time/60, time/60/60 ))

        [~,fn] = fileparts(signal_fn);
        display(sprintf('%s', fn))
        for nc = 1:Nc
            for nk = 1:Nk
                for nR = 1:NR

                    display(sprintf('qc = %g, d = %.1f, bmax = %g, D = %g ±[%g], V = %g ±[%g], k = %g ±[%g]', ...
                        seq.qc_array(nc), 2e6*R_array(nR), b_max, ...
                        fitP_array(nc,nk,nR,1), fitP_std_array(nc,nk,nR,1), ...
                        fitP_array(nc,nk,nR,2), fitP_std_array(nc,nk,nR,2), ...
                        fitP_array(nc,nk,nR,3), fitP_std_array(nc,nk,nR,3) ))
                end
            end
        end

        fitP(n_protocol,:,:,:,:) = fitP_array;
        fitP_std(n_protocol,:,:,:,:) = fitP_std_array;
        qc(n_protocol,:) = seq.qc_array;
        R(n_protocol,:) = R_array;


    end

    fitP_data.protocol = protocol;
    fitP_data.fitP = fitP;
    fitP_data.fitP_std = fitP_std;
    fitP_data.qc = qc;
    fitP_data.R = R;

    save(data_path,'fitP_data')
else

    load(data_path,'fitP_data')
end
    protocol = fitP_data.protocol;
    fitP = fitP_data.fitP;
    fitP_std = fitP_data.fitP_std;
    qc = fitP_data.qc;
    R = fitP_data.R;
    [Np, Nc, Nk, NR, ~] = size(fitP);

    
    %nc = 1:Nc; %1:Nc;
    nk = 1;
    nR = 1;
   
    x = 1:Np;
    
    figure(1),clf
    hold on 
    c = 0;
    for nc = 1:Nc
        Y = squeeze(fitP(:,nc,nk,nR,3));
        eY = squeeze(fitP_std(:,nc,nk,nR,3));
        
        errorbar(x+c, Y, eY,'.','LineWidth',3,'MarkerSize',20)
        c = c+.1;
    end
    xlim([0.5 Np+1])
    set(gca,'LineWidth',3,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',14)





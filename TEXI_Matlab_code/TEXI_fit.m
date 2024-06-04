% Fit signals from simulations.

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

root_dir = fileparts(pwd); % folder containing simulation data


fn_prepend = '';
fn_append = 'clin7';
fn_substrate_tag = '07';

simulation_data_folder = fullfile('data_TEXI','simulation');
seq_folder = 'seq_TEXI';
results_folder = fullfile('data_TEXI','results',['sim' fn_substrate_tag]);
%results_folder = fullfile('results',['sim' fn_substrate_tag]);

if (1)
    % VT only
    sim_fn_array = {...
        'signals_spheres_d1_regular_rho50',...
        'signals_spheres_d3_regular_rho50',...
        'signals_spheres_d5_regular_rho50',...
        'signals_spheres_d7_regular_rho50',...
        'signals_spheres_d9_regular_rho50',...
        'signals_spheres_d11_regular_rho50'...
        };

    % sim_fn_array = sim_fn_array(4);
    select_seq_names = {'VT_preclin7_dz1000_tm250_b2600'};


    sim_fn_array = {...
        'signals_spheres_d4_gamma_dist',...
        'signals_spheres_d8_gamma_dist'...
        };
    select_seq_names = {...
        'FEXI_preclin7_dz500_tm400_bf2000_b1300', 'VT_preclin7_dz1000_tm250_b2600'};
    % select_seq_names = select_seq_names(1);

    % ----  reduced b range  -----

    select_seq_names = {...
        'VT_preclin7_dz1000_tm250_b2600_reduceBrange'};

    sim_fn_array = {...
        'signals_spheres_d1_regular_rho50_reduceBrange',...
        'signals_spheres_d3_regular_rho50_reduceBrange',...
        'signals_spheres_d5_regular_rho50_reduceBrange',...
        'signals_spheres_d7_regular_rho50_reduceBrange',...
        'signals_spheres_d9_regular_rho50_reduceBrange',...
        'signals_spheres_d11_regular_rho50_reduceBrange',...
        'signals_spheres_d4_gamma_dist_reduceBrange',...
        'signals_spheres_d8_gamma_dist_reduceBrange'};
end

if (1)
    % only FEXI
    select_seq_names = {'FEXI_preclin7_dz500_tm400_bf2000_b1300'};
    sim_fn_array = {...
        'signals_spheres_d1_regular_rho50', ...
        'signals_spheres_d3_regular_rho50', ...
        'signals_spheres_d5_regular_rho50', ...
        'signals_spheres_d7_regular_rho50', ...
        'signals_spheres_d9_regular_rho50', ...
        'signals_spheres_d11_regular_rho50'};
    sim_fn_array = strcat(sim_fn_array,'_20231024');

end

if (1)
    % cylinders, FEXI + VT
    %select_seq_names = {'FEXI_preclin7_dz500_tm400_bf2000_b1300', 'VT_preclin7_dz1000_tm250_b2600'};
    select_seq_names = {'FEXI_preclin7_dz500_tm400_bf2000_b1300'};
    sim_fn_array = {...
        'signals_cylinders_d1_regular', ...
        'signals_cylinders_d3_regular', ...
        'signals_cylinders_d5_regular', ...
        'signals_cylinders_d7_regular', ...
        'signals_cylinders_d9_regular', ...
        'signals_cylinders_d11_regular'};
    sim_fn_array = strcat(sim_fn_array,'_rho50_20231024');
end

if (0)
    % cylinders, VT, reduced b
    select_seq_names = {'VT_preclin7_dz1000_tm250_b2600'};
    sim_fn_array = {...
        'signals_cylinders_d1_regular', ...
        'signals_cylinders_d3_regular', ...
        'signals_cylinders_d5_regular', ...
        'signals_cylinders_d7_regular', ...
        'signals_cylinders_d9_regular', ...
        'signals_cylinders_d11_regular'};
    sim_fn_array = strcat(sim_fn_array,'_rho50_20231024_reduceBrange');
end



% use Inf for no limit, only applies to Ning fit of non-FEXI data
b_max_array = 50 + [1450   1833    2217    2600];
%b_max_array = 50 + [1450   1833];
%b_limit_array = Inf;

b_min = 0*300; % test using only  b-values > b_min for FEXI detection (b2)

% 1-FEXI, 2-Ning_h, 3-Ning_h_norm, 4-Ning_Gamma, 5-Ning_g,
% 6-Ning_h3_norm, 7-Ning_g_norm, 8-Ning_h3b_norm, 9-Ning_h_intra_norm

fit_model_for_FEXI_data_array = 3; %[1 3]; %[1 3 7];
fit_model_for_nonFEXI_data_array = 3; %[3 7];

dt = 1e-5; % waveform resolution used in simulations

do_plot_loop = 0; % for debugging

do_save_fit_sig = 1;

n_inits = 10; %100; %10; % number of inital guesses


% ensure AXR model is not used for non FEXI data
if sum(ismember(fit_model_for_nonFEXI_data_array, 1)) > 0 == 1
    display('ERROR ... fit_model_for_nonFEXI_data = 1')
    return
end



for n_sim = 1:numel(sim_fn_array)
    sim_fn = sim_fn_array(n_sim);

    % load simulation
    if numel(sim_fn) > 1 sim_fn = sim_fn(1); end
    signal_fn = fullfile(root_dir, simulation_data_folder, fn_substrate_tag, char(sim_fn));

    % load simulated signals
    load(signal_fn)
    signal_seq_names = fields(signals);

    clear seq_files

    % find sequences from seq folder
    for n = 1:numel(signal_seq_names)
        seq_files(n) = dir(fullfile(root_dir, seq_folder,[signal_seq_names{n} '.mat']));
    end


    % clean up empty fields
    seq_files = seq_files(find(cellfun(@ischar,{seq_files.name})));
    seq_names = {seq_files.name};

    if ~isempty(select_seq_names)
        select_signal_seq_names = signal_seq_names(find(contains(signal_seq_names,select_seq_names)));
    end

    if isempty(fileparts(fn_substrate_tag))
        str_substrate_tag = num2str(fn_substrate_tag);
    else
        str_substrate_tag = num2str(fileparts(fn_substrate_tag));
    end


    [~, seq_fn] = fileparts(signal_fn);
    seq_fn = extractAfter(seq_fn,'signals_');
    fitName = ['sim' str_substrate_tag '_' seq_fn '_fit'];


    for n_seq = 1:numel(select_signal_seq_names)

        seq_name = select_signal_seq_names{n_seq};
        isFEXI = contains(seq_name,{'FEXI'});

        %load sequence
        file = seq_files(find(contains(seq_names, seq_name)));
        seq_fn = file.name;
        load(fullfile(file.folder, seq_fn))
        seq_fn = extractBefore(seq_fn,'.mat');


        % make arrays for fitting
        if isFEXI % b < bmax for non-FEXI sequences
            b_max_array_tmp = Inf;
            fit_model_array = fit_model_for_FEXI_data_array;
        else
            b_max_array_tmp = b_max_array;
            fit_model_array = fit_model_for_nonFEXI_data_array;
        end

        % make arrays from simulation output structure
        [sig_array, d_array, k_array, NR, Nk] = sim2array(signals, seq_name);

        R_array = 0.5*1e-6 * d_array;

        Nc = seq.nc;


        for n_fit = 1:length(fit_model_array)
            fit_model = fit_model_array(n_fit);

            for n_b_max = 1:numel(b_max_array_tmp)
                b_max = b_max_array_tmp(n_b_max);


                % ----------  fitting  ----------------------------

                % 1-FEXI, 2-Ning_h, 3-Ning_h_norm, 4-Ning_Gamma, 5-Ning_g,
                % 6-Ning_h3_norm, 7-Ning_g_norm, 8-Ning_h3b_norm, 9-Ning_h_intra_norm

                if fit_model == 1
                    fitP_array = zeros(Nc, Nk, NR, 4 + length(seq.tm_array)); % + non-filtered
                elseif fit_model == 2 || fit_model == 4 || fit_model == 5
                    fitP_array = zeros(Nc, Nk, NR, 3 + length(seq.tm_array));
                elseif fit_model == 3 || fit_model == 7
                    fitP_array = zeros(Nc, Nk, NR, 3);
                elseif fit_model == 6 || fit_model == 8 || fit_model == 9
                    fitP_array = zeros(Nc, Nk, NR, 4);
                else
                    display('ERROR ... check fitting model')
                    return
                end

                fit_sig_struct = {};
                tic

                parfor nc = 1:Nc
%                 for nc = Nc
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
                                p = Ning_h_norm_data2fit(sig', aux, opt);

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

                            fit_sig_struct(nc,nk,nR).sFIT = sFIT(:);
                            fit_sig_struct(nc,nk,nR).sub_ind = sub_ind;
                            % display(sprintf('nc/Nc = %.2f, nk/Nk = %.2f, nR/NR = %.2f',nc/Nc, nk/Nk,nR/NR))
                        end


                    end

                end
                toc


                % optionally save fit signal
                if do_save_fit_sig
                    fit_sig_array = 0*sig_array;
                    for nc = 1:Nc
                        for nk = 1:Nk
                            for nR = 1:NR
                                sub_ind = fit_sig_struct(nc,nk,nR).sub_ind;
                                sFIT = fit_sig_struct(nc,nk,nR).sFIT;
                                fit_sig_array(nk,nR,sub_ind) = sFIT;

                                sig = squeeze(sig_array(nk,nR,sub_ind));
                                sum_square_res(nc,nk,nR) = sum((sig-sFIT).^2);

                            end
                        end
                    end

                    res.sum_square_res = sum_square_res;
                    res.sig_array = sig_array;
                    res.fit_sig_array = fit_sig_array;
                end


                res.in.Nc = Nc;
                res.in.Nk = Nk;
                res.in.NR = NR;
                res.in.gc_array = seq.gc_array;
                res.in.qc_array = seq.qc_array;
                res.in.k_array = k_array;
                res.in.R_array = R_array;

                res.in.dt = dt;

                res.fitP_array = fitP_array;

                res_fn = [seq_fn '_' fitName num2str(fit_model)];

                if isFEXI
                    if b_min > 0
                        res_fn = [res_fn '_bmin' num2str(ceil(b_min))];
                    end
                else
                    if b_min > 0
                        res_fn = [res_fn '_bmin' num2str(ceil(b_min))];
                    else
                        res_fn = [res_fn '_b' num2str(ceil(max(seq.b_array(seq.b_array < b_max*1e6))*1e-6))];
                    end
                end
               
                res_fn = [res_fn '_100'];
                res.in.fn = res_fn;

                path = fullfile(root_dir, results_folder);
                warning('off', 'MATLAB:MKDIR:DirectoryExists');
                mkdir(path)
                save(fullfile(path, res_fn),'res')
                display(sprintf('%s',res_fn))

            end

        end
    end

end





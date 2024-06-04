
% combine plots with different substrates and crushers etc.
% assumes each result is for a single substrate (size)
% plotting only min and/or max crusher data
% only k vs k

% see the fitting script for more info regarding subsampling of simulation data

clear all
close all

restoredefaultpath

addpath(genpath(fullfile(pwd,'functions')));

warning('off','all')


root_dir = fileparts(fileparts(pwd));

root_path = fullfile(root_dir, 'data_TEXI', 'results','sim07');


results_case = 1; % FEXI (d3, d5), VT (d1-d11) on regular and gamma distributed spheres
fig_case = 1; % 1-8 regular spheres, 11-14 gamma distributed spheres

results_case = 2; % reducing b vange for VT on regular and gamma distributed spheres
fig_case = 1; % 1-6 regular spheres, 11-16 gamma distributed spheres

results_case = 3; % AXR fit of FEXI & Ning fit of FEXI - regular spheres
fig_case = 3; % 1-3

results_case = 4; % AXR fit of FEXI & Ning fit of FEXI - regular cylinders
fig_case = 3; % 1-3

results_case = 5; % Ning on VT + reducing b vange for VT - regular cylinders
fig_case = 1; % 1-6

results_case = 6; % Ning fit of FEXI - regular spheres (more inits)
fig_case = 2; % 1-3

results_case = 7; % Ning fit of FEXI - regular cylinders (more inits)
fig_case = 2; % 1-3

switch results_case
    case 1 % FEXI (d3, d5), VT (d1-d11) on regular and gamma distributed spheres

        FEXI_result_names_list = {...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_rho50_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_rho50_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_rho50_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_rho50_fit3'...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_rho50_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_rho50_fit7',...
            };

        VT_result_names_list = {...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d1_regular_rho50_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d1_regular_rho50_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d1_regular_rho50_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d1_regular_rho50_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d3_regular_rho50_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d3_regular_rho50_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d3_regular_rho50_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d3_regular_rho50_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d5_regular_rho50_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d5_regular_rho50_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d5_regular_rho50_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d5_regular_rho50_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d7_regular_rho50_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d7_regular_rho50_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d7_regular_rho50_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d7_regular_rho50_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d9_regular_rho50_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d9_regular_rho50_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d9_regular_rho50_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d9_regular_rho50_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d11_regular_rho50_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d11_regular_rho50_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d11_regular_rho50_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d11_regular_rho50_fit3_b2600',...
            };


        VT_result_names_list_dist = {...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit7_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit7_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit7_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d4_gamma_dist_fit7_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit7_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit7_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit7_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_spheres_d8_gamma_dist_fit7_b2600',...
            }

        switch fig_case
            case 1 % figure 1A: AXR on FEXI (var size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(0*1+[1 3]);
            case 2 % figure 1B: Ning on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(0*2+[2]);
            case 3 % figure 1C: Ning_g_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list([5]);
            case 4 % figure 2A: Ning on VT (fix size, var slice, max b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 4]);
            case 5 % figure 2B: Ning on VT (fix size, var slice, min b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 1]);
            case 6 % figure 3A: Ning on VT (var size, single slice, max b)
                fig_opt = 4;
                result_names = VT_result_names_list([4:4:end]);
            case 7 % figure 3B: Ning on VT (var size, single slice, min b)
                fig_opt = 4;
                result_names = VT_result_names_list([1:4:end]);
            case 8 % figure 4: Ning on VT (single size, single slice, var b-range)
                fig_opt = 3;
                result_names = VT_result_names_list(0*4+[1:4]);

                % ---- gamma_dist
            case 11 % figure 4A: Ning on VT (fix size distribution, var slice, max b)
                fig_opt = 3;
                result_names = VT_result_names_list_dist([8*0 + 4]);
            case 12 % figure 4B: Ning on VT (fix size distribution, var slice, min b)
                fig_opt = 3;
                result_names = VT_result_names_list_dist([8*0 + 1]);
            case 13 % figure 5A: Ning on VT (size distributions, var slice, max b)
                fig_opt = 4;
                result_names = VT_result_names_list_dist([[0 1]*8 + 4]);
            case 14 % figure 5B: Ning on VT (size distributions, var slice, min b)
                fig_opt = 4;
                result_names = VT_result_names_list_dist([[0 1]*8 + 1]);
        end






        % reducing b vange for VT
    case 2 % reducing b vange for VT on regular and gamma distributed spheres

        VT_result_names_list = {...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d1_regular_rho50_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d1_regular_rho50_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d1_regular_rho50_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d1_regular_rho50_reduceBrange_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d3_regular_rho50_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d3_regular_rho50_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d3_regular_rho50_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d3_regular_rho50_reduceBrange_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d5_regular_rho50_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d5_regular_rho50_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d5_regular_rho50_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d5_regular_rho50_reduceBrange_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d7_regular_rho50_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d7_regular_rho50_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d7_regular_rho50_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d7_regular_rho50_reduceBrange_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d9_regular_rho50_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d9_regular_rho50_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d9_regular_rho50_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d9_regular_rho50_reduceBrange_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d11_regular_rho50_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d11_regular_rho50_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d11_regular_rho50_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d11_regular_rho50_reduceBrange_fit3_b2600',...
            };

        VT_result_names_list_dist = {...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d4_gamma_dist_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d4_gamma_dist_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d4_gamma_dist_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d4_gamma_dist_reduceBrange_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d8_gamma_dist_reduceBrange_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d8_gamma_dist_reduceBrange_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d8_gamma_dist_reduceBrange_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_reduceBrange_sim07_spheres_d8_gamma_dist_reduceBrange_fit3_b2600'};



        switch fig_case
            case 1 % figure 2A: Ning on VT (fix size, var slice, max b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 4]);
                opt.append_name = '_b2600'; % add custom string for saving figures
            case 2 % figure 2B: Ning on VT (fix size, var slice, low b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 2]);
                opt.append_name = '_b1834'; % add custom string for saving figures
            case 3% figure 2B: Ning on VT (fix size, var slice, low b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 1]);
                opt.append_name = '_b1450'; % add custom string for saving figures
            case 4 % figure 3A: Ning on VT (var size, single slice, max b)
                fig_opt = 4;
                result_names = VT_result_names_list([4:4:end]);
                opt.append_name = '_b2600'; % add custom string for saving figures
            case 5 % figure 3B: Ning on VT (var size, single slice, low b)
                fig_opt = 4;
                result_names = VT_result_names_list([2:4:end]);
                opt.append_name = '_b1834'; % add custom string for saving figures
            case 6 % figure 3B: Ning on VT (var size, single slice, low b)
                fig_opt = 4;
                result_names = VT_result_names_list([1:4:end]);
                opt.append_name = '_b1450'; % add custom string for saving figures


                % ---- gamma_dist

            case 11 % figure 4A: Ning on VT (fix size distribution, var slice, max b)
                fig_opt = 3;
                result_names = VT_result_names_list_dist([4*0 + 4]);
                opt.append_name = '_b2600'; % add custom string for saving figures
            case 12 % figure 4B: Ning on VT (fix size distribution, var slice, min b)
                fig_opt = 3;
                result_names = VT_result_names_list_dist([4*0 + 2]);
                opt.append_name = '_b1834'; % add custom string for saving figures
            case 13 % figure 4B: Ning on VT (fix size distribution, var slice, min b)
                fig_opt = 3;
                result_names = VT_result_names_list_dist([4*0 + 1]);
                opt.append_name = '_b1450'; % add custom string for saving figures
            case 14 % figure 5A: Ning on VT (size distributions, var slice, max b)
                fig_opt = 4;
                result_names = VT_result_names_list_dist([[0 1]*4 + 4]);
                opt.append_name = '_b2600'; % add custom string for saving figures
            case 15 % figure 5B: Ning on VT (size distributions, var slice, low b)
                fig_opt = 4;
                result_names = VT_result_names_list_dist([[0 1]*4 + 2]);
                opt.append_name = '_b1834'; % add custom string for saving figures
            case 16 % figure 5B: Ning on VT (size distributions, var slice, low b)
                fig_opt = 4;
                result_names = VT_result_names_list_dist([[0 1]*4 + 1]);
                opt.append_name = '_b1450'; % add custom string for saving figures
        end



    case 3 % AXR fit of FEXI & Ning fit of FEXI - regular spheres

        FEXI_result_names_list = {...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d1_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d1_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d1_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d7_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d7_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d7_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d9_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d9_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d9_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d11_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d11_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d11_regular_fit7',...
            };


        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit1','_rho50_20231024_fit1');
        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit3','_rho50_20231024_fit3');
        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit7','_rho50_20231024_fit7');

        switch fig_case
            case 1 % figure 1A: AXR on FEXI (var size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(1+[1 2 3]*3);
                opt.append_name = '_AXR_d3_d5_d7';
            case 2% figure 1B: Ning_h_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(2+[1 2 3]*3);
                opt.append_name = '_Ning_d3_d5_d7';

            case 3 % figure 1C: Ning_g_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(3+[1 2 3]*3);
                opt.append_name = '_Ning_g_d3_d5_d7';
        end


    case 4 % AXR fit of FEXI & Ning fit of FEXI - regular cylinders

        FEXI_result_names_list = {...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d1_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d1_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d1_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d3_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d3_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d3_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d5_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d5_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d5_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d7_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d7_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d7_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d9_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d9_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d9_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d11_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d11_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d11_regular_fit7',...
            };


        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit1','_rho50_20231024_fit1');
        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit3','_rho50_20231024_fit3');
        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit7','_rho50_20231024_fit7');


        switch fig_case
            case 1 % figure 1A: AXR on FEXI (var size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(1+[1 2 3]*3);
                opt.append_name = '_AXR_d3_d5_d7';
            case 2 % figure 1B: Ning_h_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(2+[1 2 3]*3);
                opt.append_name = '_Ning_d3_d5_d7';

            case 3 % figure 1C: Ning_g_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(3+[1 2 3]*3);
                opt.append_name = '_Ning_g_d3_d5_d7';
        end

    case 5 % Ning on VT + reducing b vange for VT - regular cylinders

        VT_result_names_list = {...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d1_regular_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d1_regular_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d1_regular_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d1_regular_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d3_regular_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d3_regular_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d3_regular_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d3_regular_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d5_regular_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d5_regular_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d5_regular_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d5_regular_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d7_regular_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d7_regular_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d7_regular_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d7_regular_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d9_regular_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d9_regular_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d9_regular_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d9_regular_fit3_b2600',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d11_regular_fit3_b1450',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d11_regular_fit3_b1834',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d11_regular_fit3_b2217',...
            'VT_preclin7_dz1000_tm250_b2600_sim07_cylinders_d11_regular_fit3_b2600',...
            };

        VT_result_names_list = strrep(VT_result_names_list,'_fit3','_rho50_20231024_fit3');

        % for data with reduced b range apply this
        % the files with ..._reduceBrange need to be manually handled!
        if (0)
            VT_result_names_list = strrep(VT_result_names_list,'_sim07_','_reduceBrange_sim07_');
            VT_result_names_list = strrep(VT_result_names_list,'_fit3_','_reduceBrange_fit3_');
        end


        switch fig_case
            case 1  % figure 2A: Ning on VT (fix size, var slice, max b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 4]);
                opt.append_name = '_b2600'; % add custom string for saving figures
            case 2 % figure 2B: Ning on VT (fix size, var slice, low b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 2]);
                opt.append_name = '_b1834'; % add custom string for saving figures

            case 3 % figure 2B: Ning on VT (fix size, var slice, low b)
                fig_opt = 3;
                result_names = VT_result_names_list([1*4 + 1]);
                opt.append_name = '_b1450'; % add custom string for saving figures
            case 4  % figure 3A: Ning on VT (var size, single slice, max b)
                fig_opt = 4;
                result_names = VT_result_names_list([4:4:end]);
                opt.append_name = '_b2600'; % add custom string for saving figures

            case 5   % figure 3B: Ning on VT (var size, single slice, low b)
                fig_opt = 4;
                result_names = VT_result_names_list([2:4:end]);
                opt.append_name = '_b1834'; % add custom string for saving figures

            case 6 % figure 3B: Ning on VT (var size, single slice, low b)
                fig_opt = 4;
                result_names = VT_result_names_list([1:4:end]);
                opt.append_name = '_b1450'; % add custom string for saving figures

        end

    case 6 % Ning fit of FEXI - regular spheres (redo fit with more inits)

        FEXI_result_names_list = {...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d1_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d1_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d1_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d3_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d5_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d7_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d7_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d7_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d9_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d9_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d9_regular_fit7',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d11_regular_fit1',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d11_regular_fit3',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_spheres_d11_regular_fit7',...
            };

        FEXI_result_names_list = FEXI_result_names_list(contains(FEXI_result_names_list,'fit3'))


        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit1','_rho50_20231024_fit1_100');
        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit3','_rho50_20231024_fit3_100');
        FEXI_result_names_list = strrep(FEXI_result_names_list,'_fit7','_rho50_20231024_fit7_100');

        switch fig_case
            case 1 % figure 1A: AXR on FEXI (var size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(1+[1 2 3]*3);
                opt.append_name = '_AXR_d3_d5_d7';
            case 2% figure 1B: Ning_h_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list([2 3 4]);
                opt.append_name = '_Ning_d3_d5_d7_100';

            case 3 % figure 1C: Ning_g_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(3+[1 2 3]*3);
                opt.append_name = '_Ning_g_d3_d5_d7';
        end

    case 7 % Ning fit of FEXI - regular cylinders (redo fit with more inits)

        FEXI_result_names_list = {...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d1_regular_rho50_20231024_fit3_100',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d3_regular_rho50_20231024_fit3_100',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d5_regular_rho50_20231024_fit3_100',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d7_regular_rho50_20231024_fit3_100',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d9_regular_rho50_20231024_fit3_100',...
            'FEXI_preclin7_dz500_tm400_bf2000_b1300_sim07_cylinders_d11_regular_rho50_20231024_fit3_100',...
            };


        switch fig_case
            case 1 % figure 1A: AXR on FEXI (var size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(1+[1 2 3]*3);
                opt.append_name = '_AXR_d3_d5_d7';
            case 2% figure 1B: Ning_h_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list([2 3 4]);
                opt.append_name = '_Ning_d3_d5_d7_100';

            case 3 % figure 1C: Ning_g_norm on FEXI (single size)
                fig_opt = 2;
                result_names = FEXI_result_names_list(3+[1 2 3]*3);
                opt.append_name = '_Ning_g_d3_d5_d7';
        end
end
% ------------------------------------------------------

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
    seq_name = extractBefore(res_name,'_sim');
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

fig_name = strcat(sim_name,'_',seq_name,'_',substrate_name, opt.append_name);

opt.title_str = '';

opt.fig_n = 1;
fh = plot_k_list(opt);

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
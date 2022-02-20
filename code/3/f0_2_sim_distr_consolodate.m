function data_out = f0_2_sim_distr_consolodate(test_number, dgp_type, n, k_delta, k_beta2, beta1_in, beta2_in, hypothesis_type, M, array_start, array_end)
    % Consolodates the sample draws for the simulation of the robust test
    % distributions
    %
    if nargin < 1
        clear all;
        clc;
        test_number = 1;
        dgp_type = 2;
        n = 200;
        k_beta2 = 1;
        k_delta = 1;
        array_start = 1;
        array_end = 1;
        M = 100;
        beta1_in = 1 .* ones(k_delta, 1); 
        beta2_in = 0 .* ones(k_beta2, 1); 
        hypothesis_type = 1;
    end
    num_batches = array_end - array_start + 1;
        

    warning('off','all');
    
    k_lambda_n = k_beta2;

%     num_batches = 100;
%     M = 100;  

%     data_submaindir = sprintf('./cluster_sub');
    data_submaindir = sprintf('../..');

    data_maindir = sprintf('%s/data/%d/dist/n%d', data_submaindir, test_number, n);
%     
%     num_batches = 50;
%     M = 100;
%     data_maindir = sprintf('./data/cv/n%d', n);
    output_maindir = sprintf('%s/data_combined/%d/dist/n%d', data_submaindir, test_number, n);

    data_dir  = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', data_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_beta2, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000));
    outputdir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', output_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_beta2, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000));
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    
    outfilename = sprintf('%s/dist_combined_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d.mat', outputdir, test_number, dgp_type, hypothesis_type, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000), n, k_lambda_n);

    test_distr_wald = cell(k_beta2,1);
    test_distr_max = cell(k_beta2,1);
    test_distr_max_t = cell(k_beta2,1);
    
    lambda_distr_full = zeros(k_beta2, M * num_batches);
    lambda_distr_pars = zeros(k_beta2, M * num_batches);
    theta_distr_full = zeros(k_beta2*2 + k_delta + 1, M * num_batches);
    theta_distr_pars = zeros(k_beta2*(k_delta+3), M * num_batches);
    
    theta_bs2_hat_full = zeros(k_beta2*2 + k_delta + 1, M * num_batches);
    theta_bs2_hat_pars = zeros(k_beta2*(k_delta+3), M * num_batches);
    lambda_bs2_hat_full = zeros(k_beta2, M * num_batches); 
    lambda_bs2_hat_pars = zeros(k_beta2, M * num_batches);
    max_test_stat_bs2 = zeros(M * num_batches, 1);
    max_t_test_stat_bs2 = zeros(M * num_batches, 1);
    wald_Cheng_test_stat_bs2 = zeros(M * num_batches, 1);

    theta_hat_full_Taylor_bs = zeros(k_beta2*2 + k_delta + 1, M * num_batches);
    theta_hat_pars_Taylor_bs = zeros(k_beta2*2*(k_delta+2), M * num_batches);
    lambda_hat_full_Taylor_bs = zeros(k_beta2*2, M * num_batches); 
    lambda_hat_pars_Taylor_bs = zeros(k_beta2*2, M * num_batches);
    Wald_Taylor_bs = zeros(M * num_batches, 1);
    Max_Taylor_bs = zeros(M * num_batches, 1);
    Max_t_Taylor_bs = zeros(M * num_batches, 1);
    
   
    for sim_number = (array_start-1):(array_end-1)
        currentfile = sprintf('%s/dist_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d.mat', data_dir, test_number, dgp_type, hypothesis_type, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000), n, k_lambda_n, M, sim_number);
        
        fprintf('%s \n', currentfile);
        load(currentfile, 'data');
        
        for k = 1:k_beta2
            if sim_number == 0
                    test_distr_wald{k} = [];
                    test_distr_max{k} = [];
                    test_distr_max_t{k} = [];
            end
            test_distr_wald{k} = [test_distr_wald{k}; data.test_distr_wald_Cheng{k}];
            test_distr_max{k} = [test_distr_max{k}; data.test_distr_max{k}];
            test_distr_max_t{k} = [test_distr_max_t{k}; data.test_distr_max_t{k}];
        end
        ind_M = ((sim_number * M + 1) : (sim_number+1)*M);
        lambda_distr_full(:,ind_M) = data.lambda_distr_full;
        lambda_distr_pars(:,ind_M) = data.lambda_distr_pars;
        theta_distr_full(:,ind_M) = data.mathfrac_z_distr_full;
        theta_distr_pars(:,ind_M) = data.mathfrac_z_distr_pars;

        theta_bs2_hat_full(:,ind_M) = data.theta_bs2_hat_full;
        theta_bs2_hat_pars(:,ind_M) = data.theta_bs2_hat_pars;
        lambda_bs2_hat_full(:,ind_M) = data.lambda_bs2_hat_full;
        lambda_bs2_hat_pars(:,ind_M) = data.lambda_bs2_hat_pars;
        max_test_stat_bs2(ind_M) = data.max_test_stat_bs2;
        max_t_test_stat_bs2(ind_M) = data.max_t_test_stat_bs2;
        wald_Cheng_test_stat_bs2(ind_M) = data.wald_Cheng_test_stat_bs2;
        
        theta_hat_full_Taylor_bs(:,ind_M)  = data.theta_hat_full_Taylor_bs;
        theta_hat_pars_Taylor_bs(:,ind_M)  = data.theta_hat_pars_Taylor_bs;
%         lambda_hat_full_Taylor_bs(:,ind_M) = data.lambda_hat_full_Taylor_bs; 
%         lambda_hat_pars_Taylor_bs(:,ind_M) = data.lambda_hat_pars_Taylor_bs;
        Wald_Taylor_bs(ind_M)  = data.Wald_Taylor_bs;
        Max_Taylor_bs(ind_M)   = data.Max_Taylor_bs;
        Max_t_Taylor_bs(ind_M) = data.Max_t_Taylor_bs;

    end
    fprintf('\n New File: %s \n', outfilename);
    save(outfilename, 'test_distr_wald', 'test_distr_max', 'test_distr_max_t', 'lambda_distr_full', 'lambda_distr_pars', 'theta_distr_full', 'theta_distr_pars', 'theta_bs2_hat_full', 'theta_bs2_hat_pars', 'lambda_bs2_hat_full', 'lambda_bs2_hat_pars', 'max_test_stat_bs2', 'max_t_test_stat_bs2', 'wald_Cheng_test_stat_bs2', 'Wald_Taylor_bs', 'Max_Taylor_bs', 'Max_t_Taylor_bs', 'theta_hat_full_Taylor_bs', 'theta_hat_pars_Taylor_bs', 'lambda_hat_full_Taylor_bs', 'lambda_hat_pars_Taylor_bs');

end
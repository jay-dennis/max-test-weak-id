function data_out = f1_2_sims_consolodate(test_number, dgp_type, n, k_delta, k_beta2, beta1_in, beta2_in, hypothesis_type, J, array_start, array_end)
    % Consolodates the sample draws for the simulation of the robust test
    %
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
        J = 100;
        beta1_in = 1 .* ones(k_delta, 1); 
        beta2_in = 0 .* ones(k_beta2, 1); 
        hypothesis_type = 1;
    end

    warning('off','all');
    
    k_lambda_n = k_beta2;

%     data_submaindir = sprintf('./cluster_sub');
    data_submaindir = sprintf('../..');

    data_maindir = sprintf('%s/data/%d/sims/n%d', data_submaindir, test_number, n);
%     data_maindir = sprintf('./data/cv/n%d', n);
    output_maindir = sprintf('%s/data_combined/%d/sims/n%d', data_submaindir, test_number, n);

    data_dir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', data_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_beta2, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000));
    outputdir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', output_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_beta2, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000));
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    
    outfilename = sprintf('%s/sims_combined_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d.mat', outputdir, test_number, dgp_type, hypothesis_type, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000), n, k_lambda_n);

    data_all = [];
    for sim_number = (array_start-1):(array_end-1)
        currentfile = sprintf('%s/sims_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d.mat', data_dir, test_number, dgp_type, hypothesis_type, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000), n, k_lambda_n, J, sim_number);
        
        fprintf('%s \n', currentfile);
        load(currentfile, 'data0');
        
        data_all = [data_all data0];
        clear data0;
    end
    data0 = data_all;
    
    fprintf('\n New File: %s \n', outfilename);
    save(outfilename, 'data0');

end
function data_out = f1_2_sims_consolodate(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, array_start, array_end)
    % Consolodates the sample draws for the simulation of the robust test
    %
    %
    if nargin < 1
        clear all;
        clc;
        dgp_type = 2;
        n = 200;
        k_lambda = 1;
        k_delta = 1;
        array_start = 1;
        array_end = 1;
        J = 100;
        beta1 = 1; 
        beta2 = 0; 
        hypothesis_type = 1;
    end
    test_number = 1;

    warning('off','all');
    
    beta_in = [beta1 beta2];
    k_lambda_n = k_lambda;

    data_maindir = sprintf('../../data/%d/sims/n%d', test_number, n);
%     data_maindir = sprintf('./data/cv/n%d', n);
    output_maindir = sprintf('../../data_combined/%d/sims/n%d', test_number, n);

    data_dir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', data_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000));
    outputdir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', output_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000));
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    
    outfilename = sprintf('%s/sims_combined_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d.mat', outputdir, test_number, dgp_type, hypothesis_type, floor(beta1*1000), floor(beta2*1000), n, k_lambda_n);

    data_all = [];
    for sim_number = (array_start-1):(array_end-1)
        currentfile = sprintf('%s/sims_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d.mat', data_dir, test_number, dgp_type, hypothesis_type, floor(beta1*1000), floor(beta2*1000), n, k_lambda_n, J, sim_number);
        
        fprintf('%s \n', currentfile);
        load(currentfile, 'data0');
        
        data_all = [data_all data0];
        clear data0;
    end
    data0 = data_all;
    
    fprintf('\n New File: %s \n', outfilename);
    save(outfilename, 'data0');

end
function data_out = f1_2_sims_cleanup(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, array_start, array_end)
    % Consolodates the sample draws for the simulation of the robust test
    %
    %
    if nargin < 1
        clear all;
        clc;
        dgp_type = 2;
        n = 200;
        k_lambda = 20;
        k_delta = 1;
        array_start = 1;
        array_end = 100;
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

    data_dir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', data_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000));

    out_maindir = sprintf('../../data/%d/sims/temp/n%d', test_number, n);
    out_dir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', out_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000));
    if exist(out_dir,'dir') == 0
        mkdir(out_dir);
    end

    for sim_number = (array_start-1):(array_end-1)
        currentfile = sprintf('%s/sims_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d.mat', data_dir, test_number, dgp_type, hypothesis_type, floor(beta1*1000), floor(beta2*1000), n, k_lambda_n, J, sim_number);
%         outfile = sprintf('%s/sims_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d.mat', out_dir, test_number, dgp_type, hypothesis_type, floor(beta1*1000), floor(beta2*1000), n, k_lambda_n, J, sim_number);
        outfile = currentfile;
        
        fprintf('\n %s \n', currentfile);
%         load(currentfile, 'data0');
        load(currentfile);
        
        for j = 1:J
            data0(j) = data0(j).clean_up();
        end
        
        fprintf('\n New File: %s \n \n', outfile); 
        save(outfile, 'data0', 'time');

        clear data0 time;
    end

end
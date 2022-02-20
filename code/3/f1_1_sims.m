function f1_1_sims(test_number, dgp_type, n, beta1_in, beta2_in, pi_in, hypothesis_type, J, seed_type)
% function f1_1_sims(dgp_type, n, k_delta, k_beta2, beta1, beta2, hypothesis_type, J, seed_type)
% Calculates tables for the robust test


if nargin < 1
    clear all; clc;
    J = 2;            %number of Simulations
    dgp_type = 2;
    n = 200;
    seed_type = 1;

    % k_lambda_vec = [1 20];
    k_beta2 = 1;
    k_delta = 1; % the first is always an intercept when this > 1. if this=1 then no intercept

    beta1_in = 1 .* ones(k_delta, 1); 
    beta2_in = 0 .* ones(k_beta2, 1); 
    pi_in = ones(k_delta+k_beta2, 1);
    hypothesis_type = 1;
    test_number = 1;
else
    k_delta = length(beta1_in);
    k_beta2 = length(beta2_in);
end


    if seed_type == 0  % all at once, random
        seed_input = 'shuffle';
        rng('shuffle');
        temp = rng;
        sim_number = temp.Seed;
        clear temp;
    elseif seed_type == 1  % all at once, reproducible
        seed_input0 = 1;
        seed_input0 = seed_input0 - 1;
        sim_number = seed_input0;
    elseif seed_type == 2  % longleaf array, reproducible
        seed_input0 = str2num(getenv('SLURM_ARRAY_TASK_ID'));
        seed_input0 = seed_input0 - 1;
        sim_number = seed_input0;
    end
    
    k_beta2_n = k_beta2;
   
tic;    
    warning('off','all');

    for j = 1:J
        if (seed_type == 1 || seed_type == 2)
            seed_input = seed_input0 * 10^ceil(log(J)/log(10)) + j;
        end
        data0(j) = class_tests_1(test_number, dgp_type, hypothesis_type, n, seed_input, beta1_in, beta2_in, pi_in, k_beta2_n);
    end

    for j = 1:J
        data0(j) = data0(j).fcn_estimation_initial();
        data0(j) = data0(j).fcn_tests_initial();
        data0(j) = data0(j).clean_up();
    end

    
    
%%
    output_main_dir = sprintf('./data/%d/sims/n%d', test_number, n);

    outputdir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', output_main_dir, test_number, dgp_type, hypothesis_type, n, k_delta, k_beta2, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000));
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    outputname=sprintf('%s/sims_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d.mat', outputdir, test_number, dgp_type, hypothesis_type, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000), n, k_beta2_n, J, sim_number);

    count = 0;
    while exist(outputname,'file') == 2
        fprintf('Error: Name Taken; Iterating Name Forward \n')
        count = count + 1;
        outputname=sprintf('%s/sims_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d_J%d_%d_%d.mat', outputdir, test_number, dgp_type, hypothesis_type, floor(beta1_in(1)*1000), floor(beta2_in(1)*1000), n, k_beta2_n, J, sim_number, count);
    end
    fprintf('%s \n', outputname);
    time = toc / 60;
    
    save(outputname, 'data0', 'time');

    
    

function f0_sim_wrapper(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, M, seed_type, S)
    if nargin < 1
        clear; clc;
        J = 1000; % 2000;
        M = 100; % 100           %number of Simulations
        dgp_type = 2;
        n = 200;
        k_delta = 1;
        k_lambda = 2;
        seed_type = 1;
        S = 10; %20; % 100
        beta1 = 1; 
        beta2 = 0;
        hypothesis_type = 1;
    end
    
    for s = 1:S
        f0_1_sim_distr(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, M, seed_type, s);
    end
    % f0_2_sim_distr_consolodate(dgp_type, n, k_delta, k_lambda, M, S);
    f1_1_sims(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, seed_type);
end
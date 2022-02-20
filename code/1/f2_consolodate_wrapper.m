function f2_consolodate_wrapper(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, array_start, array_end)
    if nargin < 1
        clear; clc;
        J = 100; 
        dgp_type = 2;
        n = 200;
        k_delta = 1;
        k_lambda = 2;
        array_start = 1;
        array_end = 1;
        beta1 = 1; 
        beta2 = 0;
        hypothesis_type = 1;
    end
    M=J;
    % b1_vec = unique([0 1 2 5 10 1*sqrt(n)]);
    b1_vec = unique([0 1*sqrt(n)]);
    beta1_vec = b1_vec / sqrt(n);
    % b2_vec = unique([0 1 1*sqrt(n)]);
    b2_vec = unique([0 1]);
    beta2_vec = b2_vec / sqrt(n);
    for beta1 = beta1_vec
      for hypothesis_type = 1:length(b2_vec)
        beta2 = beta2_vec(hypothesis_type);
        f0_2_sim_distr_consolodate(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, M, array_start, array_end);
        f1_2_sims_consolodate(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, array_start, array_end);
      end
    end
end
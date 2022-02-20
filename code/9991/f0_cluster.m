function f0_cluster(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, seed_type)
    if nargin < 1
        clear; clc;
        J = 10; 
        dgp_type = 2;
        n = 200;
        k_delta = 1;
        k_lambda = 3;
        seed_type = 1;
        %         beta1 = 1; 
        %         beta2 = 0;
        %         hypothesis_type = 1;
    end
    M=J;
    b1_vec = unique([0 1 2 5 10 1*sqrt(n)]);
    %     b1_vec = unique([0 1 2 5 10]);
    beta1_vec = b1_vec / sqrt(n);
    b2_vec = unique([0 1 1*sqrt(n)]);
    beta2_vec = b2_vec / sqrt(n);
    for beta1 = beta1_vec
      for hypothesis_type = 2:length(b2_vec)
        beta2 = beta2_vec(hypothesis_type);
        f0_1_sim_distr(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, M, seed_type);
        f1_1_sims(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, seed_type);
      end
    end
end
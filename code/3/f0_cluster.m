function f0_cluster(test_number, dgp_type, n, k_delta, k_beta2, J, seed_type)
    if nargin < 1
        clear; clc;
        J = 500; 
        test_number = 1;
        dgp_type = 2;
        n = 200;
        k_delta = 1;
        k_beta2 = 2;
        seed_type = 1;
        %         beta1 = 1; 
        %         beta2 = 0;
        %         hypothesis_type = 1;j
    end
    M=J;
    %
    if test_number == 1
        hyp_vec = [1 2 3]; % testing beta2 here
        % b1_vec = unique([0 1 2 5 10 1*sqrt(n)]);
        b1_vec = unique([0 1*sqrt(n)]);
        b1_vec = repmat(b1_vec, k_delta, 1);
        beta1_vec = b1_vec ./ sqrt(n);
        b2_vec = unique([0 1 1*sqrt(n)]);
        b2_vec = repmat(b2_vec, k_beta2, 1);
        beta2_vec = b2_vec ./ sqrt(n);
        pi_vec = 0 .* ones(k_delta+k_beta2,1);
        pi_in = pi_vec;
    end
    if test_number == 2
        hyp_vec = [1 1 1]; % bc testing pi here under different id strengths
        b1_vec = unique([0 1*sqrt(n)]);
        b1_vec = repmat(b1_vec, k_delta, 1);
        beta1_vec = b1_vec ./ sqrt(n);
        b2_vec = unique([0 1 1*sqrt(n)]);
        b2_vec = repmat(b2_vec, k_beta2, 1);
        beta2_vec = b2_vec ./ sqrt(n);
        pi_vec = 0 .* ones(k_delta+k_beta2,1);
        pi_in = pi_vec;
    end
    if test_number == 3
        hyp_vec = [2 2 2]; % bc testing pi here under different id strengths
        b1_vec = unique([0 1*sqrt(n)]);
        b1_vec = repmat(b1_vec, k_delta, 1);
        beta1_vec = b1_vec ./ sqrt(n);
        b2_vec = unique([0 1 1*sqrt(n)]);
        b2_vec = repmat(b2_vec, k_beta2, 1);
        beta2_vec = b2_vec ./ sqrt(n);
        pi_vec = (1/sqrt(n)) .* ones(k_delta+k_beta2,1);
        pi_in = pi_vec;
    end
    if test_number == 4
        hyp_vec = [3 3 3]; % bc testing pi here under different id strengths
        b1_vec = unique([0 1*sqrt(n)]);
        b1_vec = repmat(b1_vec, k_delta, 1);
        beta1_vec = b1_vec ./ sqrt(n);
        b2_vec = unique([0 1 1*sqrt(n)]);
        b2_vec = repmat(b2_vec, k_beta2, 1);
        beta2_vec = b2_vec ./ sqrt(n);
        pi_vec = (1) .* ones(k_delta+k_beta2,1);
        pi_in = pi_vec;
    end
    if test_number == 5
        hyp_vec = 1; % bc testing pi here where pi controls id strength of beta2
        b1_vec = unique([1*sqrt(n)]);
        b1_vec = repmat(b1_vec, k_delta, 1);
        beta1_vec = b1_vec ./ sqrt(n);
        b2_vec = unique([1*sqrt(n)]);
        b2_vec = repmat(b2_vec, k_beta2, 1);
        beta2_vec = b2_vec ./ sqrt(n);
        pi_vec = 0 .* ones(k_delta+k_beta2,1);
        pi_in = pi_vec;
    end
    if test_number == 6
        hyp_vec = 1; % bc testing pi here where pi controls id strength of beta2
        b1_vec = unique([1*sqrt(n)]);
        b1_vec = repmat(b1_vec, k_delta, 1);
        beta1_vec = b1_vec ./ sqrt(n);
        b2_vec = unique([1*sqrt(n)]);
        b2_vec = repmat(b2_vec, k_beta2, 1);
        beta2_vec = b2_vec ./ sqrt(n);
        pi_vec = (1/sqrt(n)) .* ones(k_delta+k_beta2,1);
        pi_in = pi_vec;
    end
    %
    for ind_beta1 = 1:size(beta1_vec,2)
        beta1_in = beta1_vec(:,ind_beta1);
        for ind_hyp = 1:length(hyp_vec)
            hypothesis_type = hyp_vec(ind_hyp);
            beta2_in = beta2_vec(:,ind_hyp);
%             f0_1_sim_distr(dgp_type, n, k_delta, k_beta2, beta1_in, beta2_in, hypothesis_type, M, seed_type);
            f0_1_sim_distr(test_number, dgp_type, n, beta1_in, beta2_in, pi_in, hypothesis_type, M, seed_type);
            f1_1_sims(test_number, dgp_type, n, beta1_in, beta2_in, pi_in, hypothesis_type, J, seed_type);
        %   f1_1_sims(dgp_type, n, k_delta, k_beta2, beta1, beta2, hypothesis_type, J, seed_type);
        end
    end
    fprintf('\n Done! \n');
end
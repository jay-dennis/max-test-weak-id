% Temp
    clear; clc;

    J = 10;            %number of Simulations
    dgp_type = 2;
    n = 200;
    seed_type = 1;
    seed_input0 = 0;
    sim_number = seed_input0;

    test_number = 1;

    k_lambda = 20;
    k_delta = 1;
    k_lambda_n = k_lambda;

    beta1 = 1; 
    beta2 = 0; % 1/sqrt(n); 
    hypothesis_type = 1;

    beta_in = [beta1 beta2];

    warning('off','all');

    for j = 1:J
        if (seed_type == 1 || seed_type == 2)
            seed_input = seed_input0 * 10^ceil(log(J)/log(10)) + j;
        end
        data0(j) = class_tests_1(dgp_type, hypothesis_type, n, seed_input, k_lambda, k_delta, beta_in, k_lambda_n);
    end

    for j = 1:J
        data0(j) = data0(j).fcn_estimation_initial();
        data0(j) = data0(j).fcn_tests_initial();
    end

    tic;
    for j = 1:J
        data0(j) = data0(j).fcn_cluster_sim_distr;
    end
    fprintf('%.2f \n \n', toc/60);
    
    
    
    obj = data0(1);
    
    
    % test full stuff
    Y = obj.Y0;
    X = [obj.x_delta; obj.x_lambda];
    theta = obj.theta_hat_full;        
    theta = obj.theta_0;        
    i=[];
    type = 0;
    type_exclude_beta = 1;
    l_K = obj.fcn_l_K_selection(theta, X, Y, type, i, type_exclude_beta);
    
    theta_0_in = obj.theta_0;
    Z = ones(1, obj.n);
    
    tau = obj.fcn_tau(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
            
    chi = obj.fcn_chi(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
    
    pi_star = obj.fcn_pi_star(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
    
    mathfrac_z = obj.fcn_mathfrac_z(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
    
    mathfrac_z_distr = obj.fcn_mathfrac_z_distr(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in);
    
    Z = ones(1, obj.n);
    type = 0;
    mathfrac_z = obj.fcn_mathfrac_z_wrapperfcn_mathfrac_z(type, Z);
    type = 1;
    mathfrac_z = obj.fcn_mathfrac_z_wrapperfcn_mathfrac_z(type, Z);
    
    M=1;
    type = 0;
    mathfrac_z_distr_full = obj.fcn_mathfrac_z_distr(type, M);
    type = 1;
    mathfrac_z_distr_pars = obj.fcn_mathfrac_z_distr(type, M);
    
    
    tic;
    obj = obj.fcn_cluster_sim_distr;
    fprintf('%.2f \n \n', toc/60);
 
    
    
    
    
    
    
    
    
    
    
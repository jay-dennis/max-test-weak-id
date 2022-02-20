function f1_0_results_baseline_rfs()
    % Make the initial tables and histograms to show comparison with
    % standard critical values
    clear; clc;
%%
% Setup Parameters
    J = 1000;            %number of Simulations
    dgp_type = 2;
    n = 200;
    seed_type = 1;
    seed_input0 = 0;
    sim_number = seed_input0;

    test_number = 1;

%     k_lambda = 1;  % Run for k_lambda = 1, 20
    k_lambda_vec = [1 20];
    k_delta = 1;

    beta1 = 1; 
    beta2 = 0; % 1/sqrt(n); 
    hypothesis_type = 1;

    beta_in = [beta1 beta2];

    alpha_vec = [0.01 0.05 0.10];

    outputdir = '../output';

    fig = 0;
    
    for k_lambda = k_lambda_vec
            k_lambda_n = k_lambda;

%%
    % Linear Model with no weak id.
%     clear; clc;
    X = randn(n, k_delta+k_lambda);
    lambda = zeros(k_lambda,1);
    lambda(1) = beta2;
    theta = [ones(k_delta, 1); lambda];

    lambda_vec = zeros(J,k_lambda);
    lambda_vec_pars = zeros(J,k_lambda);
    wald_test = zeros(J,1);
    max_test = zeros(J,1);
    max_temp = zeros(J,1);
    for j = 1:J
        if (seed_type == 1 || seed_type == 2)
            seed_input = seed_input0 * 10^ceil(log(J)/log(10)) + j;
        end
        rng(seed_input);
        
        e = randn(n,1);
        y = X * theta + e;

        % For Wald Test
        theta_hat = X \ y;
        e_hat = y - X * theta_hat;
        lambda_hat = theta_hat(k_delta+1:end);
        lambda_vec(j,:) = lambda_hat;
        inner = X .* repmat(e_hat, 1, k_lambda+k_delta)  / sqrt(n);
        outer = inv(X' * X / n);
        V = outer * (inner' * inner) * outer';
        V_lambda = V(k_delta+1:end, k_delta+1:end);
        wald_test(j) = (n-k_lambda) * lambda_hat' * inv( V_lambda ) * lambda_hat;

        % For Max Test
        lambda_hat_pars = zeros(k_lambda,1);
        Vi = zeros((k_delta+1)*k_lambda,n);
        % se_lambda = zeros(k,1);
        for i = 1:k_lambda
            X_temp = [X(:,1:k_delta) X(:,k_delta+i)];
            theta_hat_pars_i = X_temp \ y;
            nu_hat = y - X_temp * theta_hat_pars_i;
            lambda_vec_pars(j,i) = theta_hat_pars_i(end);
            lambda_hat_pars(i) = theta_hat_pars_i(end);
            inner = X_temp .* repmat(nu_hat, 1, k_delta+1) / sqrt(n);
            outer = inv(X_temp' * X_temp / n);
            ind = [(k_delta+1)*i-k_delta : (k_delta+1)*i];
            Vi(ind,:) =  outer * inner';
        end
        V_pars = Vi * Vi';
        % se_lambda(i) = sqrt(V(end,end));
        ind = ((k_delta+1):(k_delta+1):(k_delta+1)*k_lambda);
        V_pars_lambda = V_pars(ind,ind);
        Wt = sqrt(diag(V_pars_lambda)).^(-1);
        max_test(j) = sqrt(n) * max(abs(lambda_hat_pars));
        max_t_test(j) = sqrt(n) * max(abs(Wt .* lambda_hat_pars));
        theta_sim = randn((k_delta+1)*k_lambda,1); 
        theta_sim = sqrtm(V_pars) * theta_sim;
        lambda_sim = theta_sim(ind);
        max_temp(j) = max(abs(lambda_sim))';
        max_t_temp(j) = max(abs(Wt .* lambda_sim))';

    end
        chi2_temp = chi2rnd(k_lambda,J,1);

    chi2_temp_sorted = sort(chi2_temp);
    wald_sorted = sort(wald_test);
    max_test_sorted = sort(max_test);
    max_t_test_sorted = sort(max_t_test);
    max_temp_sorted = sort(max_temp);
    max_t_temp_sorted = sort(max_t_temp);
    % [max_test_sorted max_temp_sorted]
    % [chi2_temp_sorted wald_Cheng_sorted]
    for i = 1:length(alpha_vec)
        a = alpha_vec(i);
        ind_a = floor(J*(1-a));
        wald_tail{i} = wald_sorted(ind_a:end);
        chi_tail{i} = chi2_temp_sorted(ind_a:end);
        max_test_tail{i} = max_test_sorted(ind_a:end);
        max_t_test_tail{i} = max_t_test_sorted(ind_a:end);
        max_norm_tail{i} = max_temp_sorted(ind_a:end);
        max_t_norm_tail{i} = max_t_temp_sorted(ind_a:end);
        rf_reg(:,i) = [sum((wald_sorted >= chi_tail{i}(1)))/J; sum((max_test_sorted >= max_norm_tail{i}(1)))/J; sum((max_t_test_sorted >= max_t_norm_tail{i}(1)))/J];
    end
    
    if k_lambda == 1
        windowsize_tests = [0 0 .75 .5];
        windowsize_lambda = [0 0 .25 .5];
    elseif k_lambda == 20
        windowsize_tests = [0 0 .75 .5];
        windowsize_lambda = [0 0 1 .5];
    end

    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_tests);
        subplot(1,2,1);
            histogram(chi2_temp, 'BinWidth', 1, 'Normalization', 'probability');
            hold on;
            histogram(wald_test, 'BinWidth', 1, 'Normalization', 'probability');
            hold off;
        subplot(1,2,2);
            histogram(max_temp, 'BinWidth', .1, 'Normalization', 'probability');
            hold on;
            histogram(max_test, 'BinWidth', .1, 'Normalization', 'probability');
            hold off;
%         subplot(1,3,3);
%         histogram(max_t_temp, 'BinWidth', .1, 'Normalization', 'probability');
%         hold on;
%         histogram(max_t_test, 'BinWidth', .1, 'Normalization', 'probability');
%         hold off;
        suptitle('Linear Model Tests');
%         temp_title = sprintf('$$ \\hat{\\rho}_{n} $$, %s', innovation_type_string);
%         title(temp_title, 'Interpreter', 'latex', 'FontSize', 22);
%         axis([0, 16, -.4, .4]);
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'linear_tests', J, n, k_lambda);
        saveas(fig,outputname,'png');

   
    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_lambda);
            for i = 1:k_lambda
                subplot(2,k_lambda,i);
                histogram(lambda_vec(:,i));
                subplot(2,k_lambda,k_lambda+i);
                histogram(lambda_vec_pars(:,i));
            end
        suptitle('Linear Model Lambdas');
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'linear_lambdas', J, n, k_lambda);
        saveas(fig,outputname,'png');

    %%
    % Models with weak identification

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

    theta_hat_full = zeros(length(data0(1).theta_0), J);
    theta_hat_pars = zeros(data0(1).k_delta+1+2, data0(1).k_lambda, J);
    lambda_hat_full = zeros(k_lambda, J);
    lambda_hat_pars = zeros(k_lambda, J);
    wald_Cheng = zeros(J,1);
    max_test = zeros(J,1);
    max_t_test = zeros(J,1);
    max_temp = zeros(J,1);

    for j = 1:J
        theta_hat_full(:,j) = data0(j).theta_hat_full;
        theta_hat_pars(:,:,j) = data0(j).theta_hat_pars;
        lambda_hat_full(:,j) = data0(j).lambda_hat_full;
        lambda_hat_pars(:,j) = data0(j).lambda_hat_pars;
        wald_Cheng(j) = data0(j).wald_Cheng_test_stat;
        max_test(j) = data0(j).max_test_stat;
        max_t_test(j) = data0(j).max_t_test_stat;
        max_temp(j) = data0(j).max_temp;
        max_t_temp(j) = data0(j).max_t_temp;
    end
        chi2_temp = chi2rnd(k_lambda,J,1);
        % max_temp = randn(k_lambda,J); max_temp = max(abs(max_temp), [], 1)';

    chi2_temp_sorted = sort(chi2_temp);
    wald_Cheng_sorted = sort(wald_Cheng);
    max_test_sorted = sort(max_test);
    max_temp_sorted = sort(max_temp);
    max_t_test_sorted = sort(max_t_test);
    max_t_temp_sorted = sort(max_t_temp);
    % [max_temp_sorted max_test_sorted]
    % [chi2_temp_sorted wald_Cheng_sorted]
    for i = 1:length(alpha_vec)
        a = alpha_vec(i);
        ind_a = floor(J*(1-a));
        wald_tail{i} = wald_Cheng_sorted(ind_a:end);
        chi_tail{i} = chi2_temp_sorted(ind_a:end);
        max_test_tail{i} = max_test_sorted(ind_a:end);
        max_norm_tail{i} = max_temp_sorted(ind_a:end);
        max_t_norm_tail{i} = max_t_temp_sorted(ind_a:end);
        rf_wi(:,i) = [sum((wald_Cheng_sorted >= chi_tail{i}(1)))/J; sum((max_test_sorted >= max_norm_tail{i}(1)))/J; sum((max_t_test_sorted >= max_t_norm_tail{i}(1)))/J];
    end

%     [alpha_vec; rf_wi]

%     figure(3);
%         subplot(2,3,1);
%         histogram(chi2_temp, 'Normalization', 'probability');
%     %     subplot(2,3,2);
%     %     histogram(wald_orig, 'Normalization', 'probability');
%         subplot(2,3,3);
%         histogram(wald_Cheng, 'Normalization', 'probability');
%         subplot(2,3,4);
%         histogram(max_temp, 'Normalization', 'probability');
%         subplot(2,3,5);
%         histogram(max_test, 'Normalization', 'probability');
%         subplot(2,3,6);
%         histogram(max_t_test, 'Normalization', 'probability');

    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_tests);
        subplot(1,2,1);
            histogram(chi2_temp, 'BinWidth', 5, 'Normalization', 'probability');
            hold on;
            histogram(wald_Cheng, 'BinWidth', 5, 'Normalization', 'probability');
            hold off;
        subplot(1,2,2);
            histogram(max_temp, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_test, 'BinWidth', .25, 'Normalization', 'probability');
            hold off;
%         subplot(1,3,3);
%         histogram(max_t_temp, 'BinWidth', .25, 'Normalization', 'probability');
%         hold on;
%         histogram(max_t_test, 'BinWidth', .25, 'Normalization', 'probability');
%         hold off;
        suptitle('Non-Linear Model with Weak Identification - Tests');
%         temp_title = sprintf('$$ \\hat{\\rho}_{n} $$, %s', innovation_type_string);
%         title(temp_title, 'Interpreter', 'latex', 'FontSize', 22);
%         axis([0, 16, -.4, .4]);
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'nonlinear_tests', J, n, k_lambda);
        saveas(fig,outputname,'png');

    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_lambda);
        for i = 1:k_lambda
            subplot(2,k_lambda,i);
            histogram(lambda_hat_full(i,:));
            subplot(2,k_lambda,k_lambda+i);
            histogram(lambda_hat_pars(i,:));
        end
        suptitle('Non-Linear Model with Weak Identification - Lambdas');
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'nonlinear_lambdas', J, n, k_lambda);
        saveas(fig,outputname,'png');

        
        
%%        
    rej_table{k_lambda} = [alpha_vec; rf_reg; rf_wi];
    
    end
%%
% Output Table

    clear table1 rownames colnames Title;
    
    col = 2;
    table1 = [rej_table{1}(2:4,col) rej_table{1}(5:7,col) rej_table{20}(2:4,col) rej_table{20}(5:7,col)];
    a = alpha_vec(col);

    table1
    
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    table_name = 'RF_initial';
    
    i=0;
    i = i+1; rownames{i} = sprintf('Wald Test');
    i = i+1; rownames{i} = sprintf('Max Test');
    i = i+1; rownames{i} = sprintf('Max t-Test');
    
    i=0;
    i = i+1; colnames{i} = sprintf('Linear Model, $k_{\\lambda,n}=%d$', k_lambda_vec(1));
    i = i+1; colnames{i} = sprintf('Non-linear Model, $k_{\\lambda,n}=%d$', k_lambda_vec(1));
    i = i+1; colnames{i} = sprintf('Linear Model, $k_{\\lambda,n}=%d$', k_lambda_vec(2));
    i = i+1; colnames{i} = sprintf('Non-linear Model, $k_{\\lambda,n}=%d$', k_lambda_vec(2));
    
    Caption=sprintf('Rejection Frequencies: $n=%d$, $k_{\\lambda,n}=%d$', n, k_lambda_n);
    Title{1}=sprintf('Rejection Frequencies, $J=%d$ $\\alpha = %.2f$', J, a);
%      Title{2}=sprintf('Experiment: %d, DGP: %s', test_number, dgp_type_string);
    Title{2}=sprintf('$n=%d$, $k_{\\lambda,n}=%d$', n, k_lambda_n);
    outputname = sprintf('./%s/%s_n%d', outputdir, table_name, n);
    tabletotex(table1, rownames, colnames, outputname, Title);
 
        
end
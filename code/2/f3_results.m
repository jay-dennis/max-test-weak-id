function f3_results(dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J)
% Calculates tables for the robust test


if nargin < 1
    clear all; clc;
    J = 10000;            %number of Simulations
    dgp_type = 2;
    n = 200;

    k_lambda_vec = [1 20];
    %     k_lambda = 1;
    k_delta = 1; % the first is always an intercept when this > 1. if this=1 then no intercept

    beta1 = 1; 
    beta2 = 0; 
    hypothesis_type = 1;
end
k_lambda_vec = [1 20];

test_number = 1;

% beta_in = [beta1 beta2];

alpha_vec = [0.01 0.05 0.10];

outputdir = './output';
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end

fig = 0;
windowsize_tests = [0 0 .75 .5];

b1_vec = unique([0 1 2 5 10 1*sqrt(n)]);
beta1_vec = b1_vec / sqrt(n);
b1_ind = 0;
for beta1 = beta1_vec
    b1_ind = b1_ind+1;
for k_lambda = k_lambda_vec
    k_lambda_n = k_lambda;

    data_maindir = sprintf('../../data_combined/%d/dist/n%d', test_number, n);
%     data_maindir = sprintf('./data_combined/dist/n%d', n);
    data_dir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', data_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000));
    datafilename = sprintf('%s/dist_combined_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d.mat', data_dir, test_number, dgp_type, hypothesis_type, floor(beta1*1000), floor(beta2*1000), n, k_lambda_n);
    load(datafilename);
    
    test_distr_wald = sort(test_distr_wald{k_lambda_n});
    test_distr_max = sort(test_distr_max{k_lambda_n}); 
    test_distr_max_t = sort(test_distr_max_t{k_lambda_n});
    
    wald_Cheng_test_stat_bs2 = sort(wald_Cheng_test_stat_bs2);
    max_test_stat_bs2 = sort(max_test_stat_bs2);
    max_t_test_stat_bs2 = sort(max_t_test_stat_bs2);
    
    wald_linearized_bs_sorted = sort(Wald_Taylor_bs);
    max_linearized_bs_sorted = sort(Max_Taylor_bs);
    max_t_linearized_bs_sorted = sort(Max_t_Taylor_bs);
    
    M = length(test_distr_wald);
    
    
    sim_maindir = sprintf('../../data_combined/%d/sims/n%d', test_number, n);
%     sim_maindir = sprintf('./data_combined/sims/n%d', n);
    sim_dir = sprintf('%s/output_%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d', sim_maindir, test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000));
    simfilename = sprintf('%s/sims_combined_%d_dgp%d_hyp%d_b1%d_b2%d_n%d_kln%d.mat', sim_dir, test_number, dgp_type, hypothesis_type, floor(beta1*1000), floor(beta2*1000), n, k_lambda_n);
    load(simfilename);

    
%%
    hypothesis_type_string = data0(1).hypothesis_type_string;
    current_stats = sprintf('%d_dgp%d_hyp%d_n%d_kd%d_kl%d_b1%d_b2%d_J%d', test_number, dgp_type, hypothesis_type, n, k_delta, k_lambda, floor(beta1*1000), floor(beta2*1000), J);

%%
    
    wald_Cheng = zeros(J,1);
    max_test = zeros(J,1);
    max_t_test = zeros(J,1);
    wald_linearized = zeros(J,1);
    max_test_linearized = zeros(J,1);
    max_t_test_linearized = zeros(J,1);
    max_dist_standard = zeros(J,1);
    max_t_dist_standard = zeros(J,1);

    
    lambda_hat_full = zeros(k_lambda,J);
    lambda_hat_pars = zeros(k_lambda,J);
    theta_hat_full = zeros(k_lambda*2 + k_delta + 1, J);
    theta_hat_pars = zeros(k_lambda*(k_delta + 3), J);
    lambda_hat_full_linearized = zeros(k_lambda*2,J);
    lambda_hat_pars_linearized = zeros(k_lambda*2,J);
    theta_hat_full_linearized = zeros(k_lambda*2 + k_delta + 1, J);
    theta_hat_pars_linearized = zeros(k_lambda*2*(k_delta + 2), J);
    
    wald_dist_standard = chi2rnd(k_lambda,J,1);
    for j = 1:J
        % test stats %
        wald_Cheng(j) = data0(j).wald_Cheng_test_stat;
        max_test(j) = data0(j).max_test_stat;
        max_t_test(j) = data0(j).max_t_test_stat;
        wald_linearized(j) = data0(j).Wald_Taylor;
        max_test_linearized(j) = data0(j).Max_Taylor;
        max_t_test_linearized(j) = data0(j).Max_t_Taylor;
        max_dist_standard(j) = data0(j).max_temp;
        max_t_dist_standard(j) = data0(j).max_t_temp;
        
        
        
        % parameter estimates %
        lambda_hat_full(:,j) = data0(j).lambda_hat_full;
        lambda_hat_pars(:,j) = data0(j).lambda_hat_pars;
        theta_hat_full(:,j) = data0(j).theta_hat_full;
        theta_hat_pars(:,j) = data0(j).theta_hat_pars(:);
        
        lambda_hat_full_linearized(:,j) = data0(j).lambda_hat_full_Taylor;
        lambda_hat_pars_linearized(:,j) = data0(j).lambda_hat_pars_Taylor;
        theta_hat_full_linearized(:,j) = data0(j).theta_hat_full_Taylor;
        theta_hat_pars_linearized(:,j) = data0(j).theta_hat_pars_Taylor(:);
    end

    wald_dist_standard = sort(wald_dist_standard);
    max_dist_standard = sort(max_dist_standard);
    max_t_dist_standard = sort(max_t_dist_standard);
    for i = 1:length(alpha_vec)
        a = alpha_vec(i);
        ind_a = floor(M*(1-a));
%         test_distr_wald_tail{i} = test_distr_wald_sorted(ind_a:end);
%         test_distr_max_tail{i} = test_distr_max_sorted(ind_a:end);
%         test_distr_max_t_tail{i} = test_distr_max_t_sorted(ind_a:end);
%         rf_wi1(:,i) = [sum((wald_Cheng >= test_distr_wald_tail{i}(1)))/J; sum((max_test >= test_distr_max_tail{i}(1)))/J; sum((max_t_test >= test_distr_max_t_tail{i}(1)))/J];
        rf_wi0(:,i) = [sum((wald_Cheng >= wald_dist_standard(ind_a)))/J; ...
                        sum((max_test >= max_dist_standard(ind_a)))/J; ...
                        sum((max_t_test >= max_t_dist_standard(ind_a)))/J];
        rf_wi1(:,i) = [sum((wald_Cheng >= test_distr_wald(ind_a)))/J; ...
                        sum((max_test >= test_distr_max(ind_a)))/J; ...
                        sum((max_t_test >= test_distr_max_t(ind_a)))/J];
        rf_wi2(:,i) = [sum((wald_Cheng >= wald_Cheng_test_stat_bs2(ind_a)))/J;  ...
                        sum((max_test >= max_test_stat_bs2(ind_a)))/J;  ...
                        sum((max_t_test >= max_t_test_stat_bs2(ind_a)))/J];
        rf_wi3(:,i) = [sum((wald_linearized >= wald_linearized_bs_sorted(ind_a)))/J;  ...
                        sum((max_test_linearized >= max_linearized_bs_sorted(ind_a)))/J;  ...
                        sum((max_t_test_linearized >= max_t_linearized_bs_sorted(ind_a)))/J;];
    end
    rej_table{b1_ind}{k_lambda} = [alpha_vec; rf_wi0; rf_wi1; rf_wi2; rf_wi3];
    rej_table{b1_ind}{k_lambda}
    









%%  
% Figures for Test Statistics
flag_figures = 0;
if flag_figures == 1
    fig = fig+1;
    figure(fig);
    rows = 4; cols = 3; pos = 0;
    set(fig, 'units','Normalized','Position',windowsize_tests);
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(wald_dist_standard, 'BinWidth', 1, 'Normalization', 'probability');
            hold on;
            histogram(wald_Cheng, 'BinWidth', 1, 'Normalization', 'probability');
            title('Wald Standard');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(max_dist_standard, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_test, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max Standard');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(max_t_dist_standard, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_t_test, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max-t Standard');
            hold off;
%
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(test_distr_wald, 'BinWidth', 1, 'Normalization', 'probability');
            hold on;
            histogram(wald_Cheng, 'BinWidth', 1, 'Normalization', 'probability');
            title('Wald BS1');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(test_distr_max, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_test, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max BS1');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(test_distr_max_t, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_t_test, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max-t BS1');
            hold off;
%
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(wald_Cheng_test_stat_bs2, 'BinWidth', 1, 'Normalization', 'probability');
            hold on;
            histogram(wald_Cheng, 'BinWidth', 1, 'Normalization', 'probability');
            title('Wald BS2');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(max_test_stat_bs2, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_test, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max BS2');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(max_t_test_stat_bs2, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_t_test, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max-t BS2');
            hold off;
            
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(wald_linearized_bs_sorted, 'BinWidth', 1, 'Normalization', 'probability');
            hold on;
            histogram(wald_linearized, 'BinWidth', 1, 'Normalization', 'probability');
            title('Wald Linearized');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(max_linearized_bs_sorted, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_test_linearized, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max Linearized');
            hold off;
        pos = pos+1; 
        subplot(rows,cols,pos);
            histogram(max_t_linearized_bs_sorted, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(max_t_test_linearized, 'BinWidth', .25, 'Normalization', 'probability');
            title('Max-t Linearized');
            hold off;
            
        tit0 = sprintf('Tests: Non-Linear Model, %s Hyp, k_{\\lambda, n}=%d, \\beta_{1}=%.2f', hypothesis_type_string, k_lambda_n, beta1);    
        suptitle(tit0);
%         temp_title = sprintf('$$ \\hat{\\rho}_{n} $$, %s', innovation_type_string);
%         title(temp_title, 'Interpreter', 'latex', 'FontSize', 22);
%         axis([0, 16, -.4, .4]);
        outputname=sprintf('./%s/%s_%s', outputdir, 'hist_tests', current_stats);
        saveas(fig,outputname,'png');

end % test stat histograms
        
        
flag_parameter_histograms = 0;
if flag_parameter_histograms == 1
    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_tests);
    for temp_k = 1:k_lambda_n
        subplot(2,k_lambda_n,temp_k);
            histogram(lambda_distr_full(temp_k,:), 'BinWidth', .25, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold on;
            histogram((sqrt(n).*lambda_hat_full(temp_k,:)), 'BinWidth', .25, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            histogram((sqrt(n).*lambda_bs2_hat_full(temp_k,:)), 'BinWidth', .25, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold off;
        subplot(2,k_lambda_n,k_lambda_n+temp_k);
            histogram(lambda_distr_pars(temp_k,:), 'BinWidth', .25, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold on;
            histogram((sqrt(n).*lambda_hat_pars(temp_k,:)), 'BinWidth', .25, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            histogram((sqrt(n).*lambda_bs2_hat_pars(temp_k,:)), 'BinWidth', .25, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold off;
    end
        suptitle('Non-Linear Model with Weak Identification - Lambdas');
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'nonlinear_lambdas_this_paper', J, n, k_lambda);
        saveas(fig,outputname,'png');

    temp_i_vec_pars = [k_delta k_delta+1 k_delta+2 k_delta+3];
    temp_i_vec_full = [k_delta k_delta+1 k_delta+k_lambda+1 k_delta+k_lambda+2];
%     % standardization depends on identification
    
    B_full = zeros(size(theta_hat_full,1), J);
    B_pars = zeros(size(theta_hat_pars,1), J);
    for j = 1:J
        [temp1 temp2] = data0(j).fcn_B(theta_hat_full(:,j),theta_hat_pars(:,j));
        B_full(:,j) = temp1;
        B_pars(:,j) = temp2;
    end
    Ns_full = sqrt(n) .* B_full;
    Ns_pars = sqrt(n) .* B_pars;
    theta_0_pars = [repmat(data0(1).delta_0,1,k_lambda); data0(1).lambda_0'; repmat(data0(1).pi_0(1),1,k_lambda); data0(1).pi_0(2:end)'];
    theta_0_pars = theta_0_pars(:);
    standardized_theta_hat_full = Ns_full .* (theta_hat_full - repmat(data0(1).theta_0,1,J));
    standardized_theta_hat_pars = Ns_pars .* (theta_hat_pars - repmat(theta_0_pars,1,J));
    
    
    B_full_bs2 = zeros(size(theta_bs2_hat_full,1), M);
    B_pars_bs2 = zeros(size(theta_bs2_hat_pars,1), M);
    for j = 1:M
        [temp1 temp2] = data0(1).fcn_B(theta_bs2_hat_full(:,j), theta_bs2_hat_pars(:,j));
        B_full_bs2(:,j) = temp1;
        B_pars_bs2(:,j) = temp2;
    end
    Ns_full_bs2 = sqrt(n) .* B_full_bs2;
    Ns_pars_bs2 = sqrt(n) .* B_pars_bs2;
    standardized_theta_bs2_hat_full = Ns_full_bs2 .* (theta_bs2_hat_full - repmat(data0(1).theta_0,1,M));
    standardized_theta_bs2_hat_pars = Ns_pars_bs2 .* (theta_bs2_hat_pars - repmat(theta_0_pars,1,M));
    
    
    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_tests);
    num = length(temp_i_vec_pars);
    for temp_i = 1:num
        subplot(num,2,2*(temp_i-1)+1); % remember to fix this 
            histogram(theta_distr_full(temp_i_vec_full(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold on;
            histogram(standardized_theta_hat_full(temp_i_vec_full(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            histogram(standardized_theta_bs2_hat_full(temp_i_vec_full(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold off;
        subplot(num,2,2*(temp_i-1)+2);
            histogram(theta_distr_pars(temp_i_vec_pars(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold on;
            histogram(standardized_theta_hat_pars(temp_i_vec_pars(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            histogram(standardized_theta_bs2_hat_pars(temp_i_vec_pars(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold off;
    end
        suptitle('Non-Linear Model with Weak Identification - Selected Theta');
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'nonlinear_thetas_this_paper', J, n, k_lambda);
        saveas(fig,outputname,'png');
    
%     % sort([standardized_theta_hat_full(end,:)' theta_distr_full(end,:)'],1)
%     temp_i = 1;
%     %histogram(theta_distr_full(temp_i_vec_full(temp_i),:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
%     histogram(sqrt(n)*theta_hat_full(temp_i,:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
%     hold on;
%     histogram(sqrt(n)*theta_bs2_hat_full(temp_i,:), 'BinWidth', .1, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
%     hold off;        
         
    
% end















% Taylor model
    temp_i_vec_full = [1 2 3 4];
    temp_i_vec_pars = [1 2 3 4 5 6];
    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_tests);
    num = length(temp_i_vec_pars);
    for temp_i = 1:num
        if temp_i <= length(temp_i_vec_full)
        subplot(num,2,2*(temp_i-1)+1); % remember to fix this 
            histogram(sqrt(n)*theta_hat_full_Taylor_bs(temp_i_vec_full(temp_i),:), 'BinWidth', .3, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold on;
            histogram(sqrt(n)*theta_hat_full_linearized(temp_i_vec_full(temp_i),:), 'BinWidth', .3, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold off;
        end
        subplot(num,2,2*(temp_i-1)+2);
            histogram(sqrt(n)*theta_hat_pars_Taylor_bs(temp_i_vec_pars(temp_i),:), 'BinWidth', .3, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold on;
            histogram(sqrt(n)*theta_hat_pars_linearized(temp_i_vec_pars(temp_i),:), 'BinWidth', .3, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
            hold off;
    end
        suptitle('Non-Linear Model with Weak Identification and Taylor- Selected Theta');
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'nonlinear_thetas_Taylor', J, n, k_lambda);
        saveas(fig,outputname,'png');


        
        
end % parameter histograms
    

end % k_lambda loop
end % beta1 loop

%%
% Output Table
flag_tex_table = 1;
if flag_tex_table == 1
    clear table1 rownames colnames Title;
    
    col = 2;
%     table1 = [rej_table{1}(2:end,col) rej_table{20}(2:end,col)];
    a = alpha_vec(col);

    len_b1 = length(b1_vec);
    table1 = [];
    for k_lambda = k_lambda_vec
    for b1_ind = 1:len_b1
        temp = rej_table{b1_ind}{k_lambda}(2:end,col);
        table1 = [table1 temp];
    end    
    end
    
    table1
    
    table_name = 'RF';
    
    i=0;
    i = i+1; rownames{i} = sprintf('Wald Test Standard');
    i = i+1; rownames{i} = sprintf('Max Test Standard');
    i = i+1; rownames{i} = sprintf('Max t-Test Standard');
    i = i+1; rownames{i} = sprintf('Wald Test BS1');
    i = i+1; rownames{i} = sprintf('Max Test BS1');
    i = i+1; rownames{i} = sprintf('Max t-Test BS1');
    i = i+1; rownames{i} = sprintf('Wald Test BS2');
    i = i+1; rownames{i} = sprintf('Max Test BS2');
    i = i+1; rownames{i} = sprintf('Max t-Test BS2');
    i = i+1; rownames{i} = sprintf('Wald Test Taylor');
    i = i+1; rownames{i} = sprintf('Max Test Taylor');
    i = i+1; rownames{i} = sprintf('Max t-Test Taylor');
    
    i=0;
    len_b1 = length(b1_vec);
    for b1_ind = 1:len_b1
        i=i+1;
        i1 = i; colnames{i1} = sprintf( '$b_{1}=%d$', floor(b1_vec(b1_ind)) );
        i2 = i+len_b1; colnames{i2} = sprintf( '$b_{1}=%d$', floor(b1_vec(b1_ind)) );
    end
    
%     Caption=sprintf('Rejection Frequencies: $n=%d$', n);
    Title{1}=sprintf('Rejection Frequencies, Experiment: %d, DGP: %d, Hyp: %s', test_number, dgp_type, hypothesis_type_string);
    Title{2}=sprintf('$n=%d$, $J=%d$, $\\alpha = %.2f$', n, J, a);
    Title{3}=sprintf('NoMulticolumn- \\multicolumn{1}{c}{} & \\multicolumn{%d}{c}{ $k_{\\lambda,n}=%d$ } & \\multicolumn{%d}{c}{ $k_{\\lambda,n}=%d$ }', len_b1, k_lambda_vec(1), len_b1, k_lambda_vec(2));
    outputname = sprintf('./%s/%s_%s', outputdir, table_name, current_stats);
    tabletotex(table1, rownames, colnames, outputname, Title);
 
end % Latex Output
    
        
        
        

  
end
%%
% sims
% clear all; clc;

obj = data0(1);
[obj.x_delta' obj.x_beta2' obj.x_pi']
obj.theta_0
obj.theta_hat_full
obj.theta_hat_pars
theta_0_pars = zeros(size(obj.theta_hat_pars,1),size(obj.theta_hat_pars,2));
theta_0_pars(1:obj.k_delta,:) = repmat(obj.theta_0(1:obj.k_delta),1,obj.k_beta2);
theta_0_pars(obj.k_delta+1,:) = obj.theta_0(obj.k_delta+1:obj.k_delta+obj.k_beta2);
theta_0_pars(obj.k_delta+2,:) = repmat(obj.theta_0(obj.k_delta+obj.k_beta2+1),1,obj.k_beta2);
theta_0_pars(obj.k_delta+3,:) = obj.theta_0(obj.k_delta+obj.k_beta2+1+1:obj.k_delta+obj.k_beta2+1+obj.k_beta2);


J = length(data0);
T_wald = zeros(J,1);
T_max = zeros(J,1);
T_max_t = zeros(J,1);
normed_theta_full = zeros(size(data0(1).theta_hat_full,1),J);
normed_theta_pars = zeros(size(data0(1).theta_hat_pars,1),size(data0(1).theta_hat_pars,2),J);
for j = 1:J
    normed_theta_full(:,j) = data0(j).Ns_full .* (data0(j).theta_hat_full - data0(j).theta_0);
    normed_theta_pars(:,:,j) = data0(j).Ns_pars .* (data0(j).theta_hat_pars - theta_0_pars);
    T_wald(j) = data0(j).wald_Cheng_test_stat;
    T_max(j) = data0(j).max_test_stat;
    T_max_t(j) = data0(j).max_t_test_stat;
end

%size(theta_full,1)
num1 = size(normed_theta_pars,1);
num2 = obj.k_beta2; %size(theta_pars,2);
num3 = 1+obj.k_beta2; %size(obj.x_pi,1);
temp1 = [1 2 4 5];
temp2 = [1 3 4 6];
temp1 = [1 2 4 5 6 7];
temp2 = [1 3 4 5 8 9];
for k = 1:num2
    %     temp = [1 1+k 1+num2+(1:num3) 1+num2+num3*k+(1:num3)]
    temp = [1 1+k 1+num2+1 1+num2+1+k]
end

% obj.Ns = [sqrt(obj.n) .* ones(obj.k_delta+obj.k_beta2,1); ones(1+obj.k_beta2,1)];
% obj.Ns = sqrt(obj.n) .* [ones(obj.k_delta+obj.k_beta2,J); theta_full(1:obj.k_delta+obj.k_beta2,:)];
% obj.Ns_lambda = obj.Ns(obj.ind_lambda_full,:);
            
for k = 1:num2
    figure(k);
    temp = [1 1+k 1+num2+1 1+num2+1+k];
    for p = 1:num1
        [2*p-1 2*p];
        subplot(num1,2,2*p-1);
        hold on;
        % correct the normalization for pi parameters; (either nothing or sqrt(n) * beta_hat)
        % correct the normalization for pi parameters; (either nothing or sqrt(n) * beta_hat)
        % correct the normalization for pi parameters; (either nothing or sqrt(n) * beta_hat)
        % correct the normalization for pi parameters; (either nothing or sqrt(n) * beta_hat)
        % correct the normalization for pi parameters; (either nothing or sqrt(n) * beta_hat)
        % probably will need to be nothing since I'm not bootstrapping the
        % distribution for every beta_hat simulation
        % histogram(obj.Ns(temp(p)) .* (theta_full(temp(p),:) - obj.theta_0(temp(p))) + obj.theta_0(temp(p)), 'BinWidth', .1, 'Normalization', 'probability');
        histogram(normed_theta_full(temp(p),:), 'BinWidth', .1, 'Normalization', 'probability');
        % histogram(theta_full(temp(p),:));
        title(sprintf('%.1f', obj.theta_0(temp(p))));
        hold off;
        subplot(num1,2,2*p);
        hold on;
        % histogram(obj.Ns(temp(p)) .* (theta_pars(p,k,:) - obj.theta_0(temp(p))) + obj.theta_0(temp(p)), 'BinWidth', .1, 'Normalization', 'probability');
        histogram(normed_theta_pars(p,k,:), 'BinWidth', .1, 'Normalization', 'probability');
        % histogram(theta_pars(p,k,:));
        title(sprintf('%.1f', obj.theta_0(temp(p))));
        hold off;
    end
end

%p=4;
%histogram(sqrt(obj.n) .* theta_full(temp(p),:), 'BinWidth', 4, 'Normalization', 'probability');
%title(sprintf('%.1f', obj.theta_0(temp(p))));


%%
% distr sim
% clear; clc;
cluster_flag = 1;

if cluster_flag == 0
    obj = data;
    obj.theta_0
    obj.theta_hat_full
    obj.theta_hat_pars

    size(obj.mathfrac_z_distr_full)
    size(obj.mathfrac_z_distr_pars)
end

num0 = obj.k_delta+obj.k_beta2+1+obj.k_beta2;
num1 = obj.k_delta+1+1+1;
num2 = obj.k_beta2;
for k = 1:num2
    temp = [1 1+k 1+num2+1 1+num2+1+k]
    temp2= [num1*(k-1)+(1:num1)]
end

for k = 1:num2
    figure(k);
    temp = [1 1+k 1+num2+1 1+num2+1+k];
    temp2= num1*(k-1)+(1:num1);
    for p = 1:num1
        [2*p-1 2*p];
        subplot(num1,2,2*p-1);
        hold on;
% %         histogram(obj.mathfrac_z_distr_full(temp(p),:) + obj.Ns(temp(p)) .* obj.theta_0(temp(p)), 'BinWidth', .1, 'Normalization', 'probability');
        if cluster_flag == 0
        histogram(obj.mathfrac_z_distr_full(temp(p),:), 'BinWidth', .1, 'Normalization', 'probability');
        elseif cluster_flag == 1
        histogram(theta_distr_full(temp(p),:), 'BinWidth', .1, 'Normalization', 'probability');
        end
        title(sprintf('%.1f', obj.theta_0(temp(p))));
        hold off;
        subplot(num1,2,2*p);
        hold on;
% %         histogram(obj.mathfrac_z_distr_pars(temp2(p),:) + obj.Ns(temp(p)) .* obj.theta_0(temp(p)), 'BinWidth', .1, 'Normalization', 'probability');
        if cluster_flag == 0
        histogram(obj.mathfrac_z_distr_pars(temp2(p),:), 'BinWidth', .1, 'Normalization', 'probability');
        elseif cluster_flag == 1
        histogram(theta_distr_pars(temp2(p),:), 'BinWidth', .1, 'Normalization', 'probability');
        end
        title(sprintf('%.1f', obj.theta_0(temp(p))));
        hold off;
    end
end




figure(51);
subplot(3,1,1);
    histogram(T_wald, 'BinWidth', .1, 'Normalization', 'probability');
    hold on;
    histogram(test_distr_wald{obj.k_beta2}, 'BinWidth', .1, 'Normalization', 'probability');
    hold off;
subplot(3,1,2);
    histogram(T_max, 'BinWidth', .1, 'Normalization', 'probability');
    hold on;
    histogram(test_distr_max{obj.k_beta2}, 'BinWidth', .1, 'Normalization', 'probability');
    hold off;
subplot(3,1,3);
    histogram(T_max_t, 'BinWidth', .1, 'Normalization', 'probability');
    hold on;
    histogram(test_distr_max_t{obj.k_beta2}, 'BinWidth', .1, 'Normalization', 'probability');
    hold off;










%%
k=2;
p=3;

temp = [1 1+k 1+num2+1 1+num2+1+k];
temp2= num1*(k-1)+(1:num1);

figure(101);
histogram(normed_theta_full(temp(p),:), 'BinWidth', .1, 'Normalization', 'probability');
hold on;
if cluster_flag == 0
    histogram(obj.mathfrac_z_distr_full(temp(p),:), 'BinWidth', .1, 'Normalization', 'probability');
elseif cluster_flag == 1
    histogram(theta_distr_full(temp(p),:), 'BinWidth', .1, 'Normalization', 'probability');
end
hold off;

figure(102);
histogram(normed_theta_pars(p,k,:), 'BinWidth', .1, 'Normalization', 'probability');
hold on;
if cluster_flag == 0
    histogram(obj.mathfrac_z_distr_pars(temp2(p),:), 'BinWidth', .1, 'Normalization', 'probability');
elseif cluster_flag == 1
    histogram(theta_distr_pars(temp2(p),:), 'BinWidth', .1, 'Normalization', 'probability');
end
hold off;

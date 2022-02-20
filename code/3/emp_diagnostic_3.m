clear all;
clc;

J = 2000; 
test_number = 2;
dgp_type = 2;
n = 500;
k_delta = 1;
k_beta2 = 3;
seed_type = 1;
hyp_vec = [1 1 1]; % bc testing pi here under different id strengths
b1_vec = unique([0 1*sqrt(n)]);
b1_vec = repmat(b1_vec, k_delta, 1);
beta1_vec = b1_vec ./ sqrt(n);
b2_vec = unique([0 1 1*sqrt(n)]);
b2_vec = repmat(b2_vec, k_beta2, 1);
beta2_vec = b2_vec ./ sqrt(n);
pi_vec = 0 .* ones(k_delta+k_beta2,1);
pi_in = pi_vec;

ind_beta1 = 1;
beta1_in = beta1_vec(:,ind_beta1);
ind_hyp = 3;
hypothesis_type = hyp_vec(ind_hyp);
beta2_in = beta2_vec(:,ind_hyp);
            
k_delta = length(beta1_in);
k_beta2 = length(beta2_in);
seed_input0 = 1;
seed_input0 = seed_input0 - 1;
sim_number = seed_input0;
k_beta2_n = k_beta2;

warning('off','all');

for j=1:J
    seed_input = seed_input0 * 10^ceil(log(J)/log(10)) + j;
    data0(j) = class_tests_1(test_number, dgp_type, hypothesis_type, n, seed_input, beta1_in, beta2_in, pi_in, k_beta2_n);
    data0(j) = data0(j).fcn_estimation_initial();
    data0(j) = data0(j).fcn_tests_initial();
end


theta_full = zeros(size(data0(1).theta_hat_full,1),J);
theta_pars = zeros(size(data0(1).theta_hat_pars,1),size(data0(1).theta_hat_pars,2),J);
for j = 1:J
    theta_full(:,j) = data0(j).theta_hat_full;
    theta_pars(:,:,j) = data0(j).theta_hat_pars;
end


obj = data0(1);
temp = [1 2 5 6];
k=1;
p=2;
figure(101);
histogram(obj.Ns(temp(p)) .*(theta_full(temp(p),:) - obj.theta_0(temp(p))), 'BinWidth', .1, 'Normalization', 'probability');
hold on;
histogram(obj.Ns(temp(p)) .* (theta_full(temp(p),:)), 'BinWidth', .1, 'Normalization', 'probability');
hold off;

mean(obj.Ns(temp(p)) .*(theta_full(temp(p),:) - obj.theta_0(temp(p))))
mean(obj.Ns(temp(p)) .* (theta_full(temp(p),:)))

figure(102);
histogram(obj.Ns(temp(p)) .* (theta_pars(p,k,:) - obj.theta_0(temp(p))), 'BinWidth', .1, 'Normalization', 'probability');
hold on;
histogram(obj.Ns(temp(p)) .* (theta_pars(p,k,:)), 'BinWidth', .1, 'Normalization', 'probability');
hold off;

mean(obj.Ns(temp(p)) .* (theta_pars(p,k,:) - obj.theta_0(temp(p))))
mean(obj.Ns(temp(p)) .* (theta_pars(p,k,:)))

corr(obj.x_beta2')






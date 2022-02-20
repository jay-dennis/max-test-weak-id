
obj = data0(1);
[obj.x_delta' obj.x_lambda' obj.x_pi']
obj.theta_0

theta_full = zeros(size(data0(1).theta_hat_full,1),J);
theta_pars = zeros(size(data0(1).theta_hat_pars,1),size(data0(1).theta_hat_pars,2),J);
for j = 1:J
    theta_full(:,j) = data0(j).theta_hat_full;
    theta_pars(:,:,j) = data0(j).theta_hat_pars;
end

%size(theta_full,1)
num1=size(theta_pars,1);
num2=size(theta_pars,2);
num3 = size(obj.x_pi,1);
temp1 = [1 2 4 5];
temp2 = [1 3 4 6];
temp1 = [1 2 4 5 6 7];
temp2 = [1 3 4 5 8 9];
for k = 1:num2
    %     temp = [1 1+k 1+num2+(1:num3) 1+num2+num3*k+(1:num3)]
    temp = [1 1+k 1+num2+1 1+num2+1+k]
end

for k = 1:num2
    figure(k);
    temp = [1 1+k 1+num2+1 1+num2+1+k];
    for p = 1:num1
        [2*p-1 2*p];
        subplot(num1,2,2*p-1);
        % histogram(sqrt(obj.n) .* (theta_full(temp(p),:) - obj.theta_0(temp(p))) + obj.theta_0(temp(p)));
        histogram(theta_full(temp(p),:));
        title(sprintf('%.1f', obj.theta_0(temp(p))));
        subplot(num1,2,2*p);
        % histogram(sqrt(obj.n) .* (theta_pars(p,k,:) - obj.theta_0(temp(p))) + obj.theta_0(temp(p)));
        histogram(theta_pars(p,k,:));
        title(sprintf('%.1f', obj.theta_0(temp(p))));
    end
end

%p=4;
%histogram(sqrt(obj.n) .* theta_full(temp(p),:), 'BinWidth', 4, 'Normalization', 'probability');
%title(sprintf('%.1f', obj.theta_0(temp(p))));




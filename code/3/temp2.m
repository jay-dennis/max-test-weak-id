clear; clc;
B = 1;
N_vec = [100 200 500 1000];
J=5000;

B_hat = zeros(J, length(N_vec), 3);
V_hat = zeros(J, length(N_vec), 3);

rng(0);
for in = 1:length(N_vec)
    n = N_vec(in);
    for j = 1:J
        X = randn(n,1);
        e0 = randn(n,1);
        e1 = sqrt(n) .* e0;
        e2 = n .* e0;
        y0 = X*B + e0;
        y1 = X*B + e1;    
        y2 = X*B + e2;    
        B_hat(j,in,1) = X \ y0;
        B_hat(j,in,2) = X \ y1;
        B_hat(j,in,3) = X \ y2;
        V_hat(j,in,1) = inv(X'*X)*X'*(e0'*e0)*X*inv(X'*X);
        V_hat(j,in,2) = inv(X'*X)*X'*(e1'*e1)*X*inv(X'*X);
        V_hat(j,in,3) = inv(X'*X)*X'*(e2'*e2)*X*inv(X'*X);
    end
end

figure(1);
for in = 1:length(N_vec)
    n = N_vec(in);
    subplot(3,length(N_vec),in);
    histogram(sqrt(n) .* (B_hat(:,in,1) -1) );
    subplot(3,length(N_vec),in + length(N_vec));
    histogram(sqrt(n) .* (B_hat(:,in,2) -1) );
    subplot(3,length(N_vec),in + (2*length(N_vec)));
    histogram(sqrt(n) .* (B_hat(:,in,3) -1) );
end

figure(2);
for in = 1:length(N_vec)
    n = N_vec(in);
    subplot(3,length(N_vec),in);
    histogram(sqrt(n) .* (V_hat(:,in,1) -1) );
    subplot(3,length(N_vec),in + length(N_vec));
    histogram(sqrt(n) .* (V_hat(:,in,2) - sqrt(n)) );
    subplot(3,length(N_vec),in + (2*length(N_vec)));
    histogram(sqrt(n) .* (V_hat(:,in,3) - n) );
end



obj = data0(1)
X0 = [obj.x_delta; obj.x_lambda; obj.x_pi];
X_in = X0;

obj.theta_0(1) = 10;
obj.theta_0(2) = 0; obj.theta_0(6) = 0;
X0 = [(-2:4/(obj.n-1):2); obj.x_lambda; ones(1,obj.n); .1.*ones(1,obj.n)]
Y = obj.fcn_Yhat(obj.theta_0, X0, 0);
dY = obj.fcn_dYhat(obj.theta_0, X0, 0);

size(X0)
size(Y)
size(dY)
dY'

figure(1);
plot(Y);
plot(dY');
plot(X0(1,:), Y, '.');

temp = randn(1,obj.n);
figure(2);
plot(temp, (1:obj.n), '.');

function out = fcn_F(obj, pi, X_pi)
    out = pi' * X_pi;
end
function out = fcn_g(obj, X2, c, gamma) %g_E2 swaps c and pi %fcn_g_E2
    out =  1 - exp( -gamma .* (X2 - c).^2 );
end
function out = fcn_dg(obj, X2, c, gamma)
    out = ( ((X2 - c).^2) .* exp( -gamma .* (X2 - c).^2 ) );
end
function Yhat = fcn_Yhat_full(obj, theta_in, X_in)
    X_d = X_in(obj.ind_k_delta, :);
    X_l = X_in(obj.ind_k_lambda_full, :);
    X_p = X_in(obj.ind_kX_pi_full, :);
    delta_in = theta_in(obj.ind_k_delta);
    lambda_in = theta_in(obj.ind_k_lambda_full);
    pi_in = theta_in(obj.ind_k_pi_full);
    len_l = obj.k_lambda;
    for i = 1:obj.k_delta
        X_d(i,:) = obj.fcn_g(X_d(i,:), 0, delta_in(i));
    end
    for i = 1:len_l
        X_l(i,:) = obj.fcn_g(X_l(i,:), 0, lambda_in(i));
    end
    X_all = [X_d; X_l];
    F = zeros(obj.k_delta+len_l, obj.n);
    for i = 1:obj.k_delta+len_l
        ind = ( (i-1)*obj.k_pi+1 : i*obj.k_pi );
        F(i,:) = obj.fcn_F(pi_in(ind), X_p);
    end
    Yhat = sum(F .* X_all,1);
end
function dYhat = fcn_dYhat_full(obj, theta_in, X_in)
    X_d = X_in(obj.ind_k_delta, :);
    X_l = X_in(obj.ind_k_lambda_full, :);
    X_p = X_in(obj.ind_kX_pi_full, :);
    delta_in = theta_in(obj.ind_k_delta);
    lambda_in = theta_in(obj.ind_k_lambda_full);
    pi_in = theta_in(obj.ind_k_pi_full);
    len_l = obj.k_lambda;
    % g
    for i = 1:obj.k_delta
        X_d(i,:) = obj.fcn_g(X_d(i,:), 0, delta_in(i));
    end
    for i = 1:len_l
        X_l(i,:) = obj.fcn_g(X_l(i,:), 0, lambda_in(i));
    end
    g = [X_d; X_l];
    % dg
    X_d = X_in(obj.ind_k_delta, :);
    X_l = X_in(obj.ind_k_lambda_full, :);
    for i = 1:obj.k_delta
        X_d(i,:) = obj.fcn_dg(X_d(i,:), 0, delta_in(i));
    end
    for i = 1:len_l
        X_l(i,:) = obj.fcn_dg(X_l(i,:), 0, lambda_in(i));
    end
    dg = [X_d; X_l];
    % F
    F = zeros(obj.k_delta+len_l, obj.n);
    for i = 1:obj.k_delta+len_l
        ind = ( (i-1)*obj.k_pi+1 : i*obj.k_pi );
        F(i,:) = obj.fcn_F(pi_in(ind), X_p);
    end
    % dF
    dF = X_p;
    %
    gdF = zeros(obj.k_pi*(obj.k_delta+len_l),obj.n);
    for i = 1:obj.k_delta+len_l
        ind = ( (i-1)*obj.k_pi+1 : i*obj.k_pi );
        gdF(ind,:) = dF .* repmat(g(i,:), obj.k_pi, 1);
    end
    %
    dYhat = [F .* dg; gdF];
end
 
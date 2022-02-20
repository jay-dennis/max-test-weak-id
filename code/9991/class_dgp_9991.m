classdef class_dgp_9991
    properties (SetAccess=public, GetAccess=public)
	 % dgp variables
     test_number
     dgp_type, dgp_type_string
     hypothesis_type, hypothesis_type_string
     g_type, g_type_string
     n, init_n
	 Y0
     e
     k_delta, k_lambda, k_pi, k_pi_per
     ind_k_delta, ind_k_lambda_full, ind_k_pi_full, ind_k_lambda_pars, ind_k_pi_pars
     ind_kX_pi_full, ind_kX_pi_pars
     x_delta, x_lambda, x_pi
     delta_0, beta_0, lambda_0, pi_0, theta_0
     b1_0, b2_0
     % estimation variables
     kappa_n
     LB_pi, UB_pi
     theta_hat_full, theta_hat_pars, lambda_hat_full, lambda_hat_pars
     % simulation variables
     cluster_sim_distr_seed
     mathfrac_z_distr_full, mathfrac_z_distr_pars, lambda_distr_full, lambda_distr_pars
     test_distr_wald_Cheng, test_distr_max, test_distr_max_t
     % BS2 variables
     cluster_bs2_distr_seed
     theta_bs2_hat_full, theta_bs2_hat_pars, lambda_bs2_hat_full, lambda_bs2_hat_pars
     max_test_stat_bs2, max_t_test_stat_bs2, wald_Cheng_test_stat_bs2       
     % general variables
     seed_input
     rng_seed
     options_fmincon, options_lsqnonlin    
    end
    methods
        function obj = class_dgp_9991(dgp_type, hypothesis_type, n, seed_input, k_lambda, k_delta, beta_in)
            obj.test_number = 9991;
            obj.dgp_type = dgp_type;
            obj.hypothesis_type = hypothesis_type;
            obj.n = n;
            obj.k_delta = k_delta; 
            obj.k_lambda = k_lambda;
            obj.k_pi_per = 1;
            obj.k_pi = (obj.k_lambda+obj.k_delta) * obj.k_pi_per;
            %
            obj.ind_k_delta = (1:obj.k_delta);
            obj.ind_k_lambda_full = (obj.k_delta+1:obj.k_delta+obj.k_lambda);
            % obj.ind_k_pi_full = (obj.k_delta+obj.k_lambda+1:obj.k_delta+obj.k_lambda+obj.k_pi*(obj.k_delta+obj.k_lambda));
            obj.ind_k_pi_full = (obj.k_delta+obj.k_lambda+1:obj.k_delta+obj.k_lambda+obj.k_pi);
            obj.ind_kX_pi_full = (obj.k_delta+obj.k_lambda+1:obj.k_delta+obj.k_lambda+obj.k_pi);
            obj.ind_k_lambda_pars = obj.k_delta+1;
            obj.ind_k_pi_pars = (obj.k_delta+1+1:obj.k_delta+1+obj.k_pi_per*(obj.k_delta+1));
            obj.ind_kX_pi_pars = (obj.k_delta+1+1:obj.k_delta+1+obj.k_pi_per*(obj.k_delta+1));
            %
            obj.seed_input = seed_input;
            obj.kappa_n = sqrt(log(obj.n)); % for ICS cvs
            %             obj.options_fmincon = optimoptions(@fmincon,'Display','iter-detailed','Diagnostics','off','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',10000,'MaxIterations',5000,'StepTolerance',1e-20,'ConstraintTolerance',1e-20,'OptimalityTolerance',1e-20); 
            %             obj.options_lsqnonlin = optimoptions(@lsqnonlin,'Display','iter-detailed','Diagnostics','off','SpecifyObjectiveGradient',true,'Algorithm','levenberg-marquardt','MaxIterations',1000,'MaxFunctionEvaluations',2400,'StepTolerance',1e-15,'FunctionTolerance',1e-15); 
            % obj.options_fmincon = optimoptions(@fmincon,'Display','off','Diagnostics','off','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',5000,'MaxIterations',2000,'StepTolerance',1e-15,'ConstraintTolerance',1e-15,'OptimalityTolerance',1e-15); 
            obj.options_fmincon = optimoptions(@fmincon,'Display','off','Diagnostics','off','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',5000,'MaxIterations',2000,'StepTolerance',1e-15,'ConstraintTolerance',1e-15,'OptimalityTolerance',1e-15); 
            obj.options_lsqnonlin = optimoptions(@lsqnonlin,'Display','off','Diagnostics','off','SpecifyObjectiveGradient',true,'Algorithm','levenberg-marquardt','MaxIterations',1000,'MaxFunctionEvaluations',2400,'StepTolerance',1e-15,'FunctionTolerance',1e-15); 
            % lsqnonlin(fun,x0,lb,ub,options)
            %
            % Now generate the data
	        obj.init_n = obj.n + 0;
            rng(0); % reproduce the X's
            %             obj.x_pi = [ones(1, obj.n); randn(obj.k_pi-1, obj.n)]; 
            obj.x_pi = randn(obj.k_pi, obj.n); 
            switch dgp_type
                case 1
                    obj.dgp_type_string = 'Case 1: iid';
                    obj.x_delta = 1*randn(obj.k_delta,obj.n);
                    obj.x_lambda = 1*randn(obj.k_lambda,obj.n);
                case 2 
                    obj.dgp_type_string = 'Case 2: blockwise independence';
                    w_delta = randn(obj.k_delta,obj.n);
                    w_beta = randn(obj.k_lambda,obj.n);
                    rho = .5;
                    num = obj.k_delta;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num)); 
                    Sig = chol(Sig);
                    obj.x_delta = Sig * w_delta;
                    num = obj.k_lambda;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num)); 
                    Sig = chol(Sig);
                    obj.x_lambda = Sig * w_beta;
                case 3 
                    obj.dgp_type_string = 'Case 3: with-in and cross-block dependence';
                    num = obj.k_delta + obj.k_lambda;
                    temp = randn(num,obj.n);
                    rho = .5;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num));
                    Sig = chol(Sig);
                    X = Sig * temp;
                    obj.x_delta = X(1:obj.k_delta,:);
                    obj.x_lambda = X(obj.k_delta+1:end,:);
            end
%             Just replace 1st variable instead            
%             obj.x_delta = [ones(1, obj.n); obj.x_delta]; % add in an intercept
%             obj.k_delta = obj.k_delta+1; % and increase k_delta to account for the intercept
            rng(seed_input); % generate the epsilons randomly
            temp = rng;
            obj.rng_seed = temp.Seed;
            clear temp;
            obj.e = randn(1,obj.init_n);
                        
            switch hypothesis_type
                case 1
                    obj.hypothesis_type_string = 'Null';
                case 2
                    obj.hypothesis_type_string = 'Local Alternative';
                case 3
                    obj.hypothesis_type_string = 'Alternative';
            end
            
            obj.lambda_0 = beta_in(obj.k_delta+1:obj.k_delta+obj.k_lambda);

            obj.beta_0 = beta_in(1:obj.k_delta);
            obj.delta_0 = beta_in(1:obj.k_delta);

            rng(0);
            % temp = repmat([1; zeros(obj.k_pi,1)],obj.k_lambda+obj.k_delta,1);
            % obj.pi_0 = temp(1:(obj.k_lambda+obj.k_delta)*obj.k_pi);
            obj.pi_0 = ones(obj.k_pi,1);
            % obj.pi_0 = [0; 1];
            
            obj.b1_0 = beta_in(1:obj.k_delta) * sqrt(obj.n);
            obj.b2_0 = obj.lambda_0 * sqrt(obj.n);
                                  
            obj.theta_0 = [obj.delta_0; obj.lambda_0; obj.pi_0];
            X0 = [obj.x_delta; obj.x_lambda; obj.x_pi];

            Ytemp = obj.fcn_Yhat(obj.theta_0, X0, 0);
            obj.Y0 = Ytemp + obj.e;
            obj.Y0 = obj.Y0((obj.init_n - obj.n + 1):end); % remove Y0 burn-in values
        end
        %
        function Yhat = fcn_Yhat(obj, theta_in, X_in, model_type)
            if model_type == 0
                Yhat = obj.fcn_Yhat_full(theta_in, X_in);
            elseif model_type == 1
                Yhat = obj.fcn_Yhat_pars(theta_in, X_in);
            end
        end
        function Yhat = fcn_Yhat_pars(obj, theta_in, X_in)
            X_d = X_in(obj.ind_k_delta,:);
            X_l = X_in(obj.ind_k_lambda_pars,:);
            X_p = X_in(obj.ind_kX_pi_pars,:);
            delta_in = theta_in(obj.ind_k_delta);
            lambda_in = theta_in(obj.ind_k_lambda_pars);
            pi_in = theta_in(obj.ind_k_pi_pars);
            len_l = 1;
            for i = 1:obj.k_delta
                X_d(i,:) = obj.fcn_g(X_d(i,:), 0, delta_in(i));
            end
            for i = 1:len_l
                X_l(i,:) = obj.fcn_g(X_l(i,:), 0, lambda_in(i));
            end
            X_all = [X_d; X_l];
            F = zeros(obj.k_delta+len_l, obj.n);
            for i = 1:obj.k_delta+len_l
                ind = ( (i-1)*obj.k_pi_per+1 : i*obj.k_pi_per );
                F(i,:) = obj.fcn_F(pi_in(ind), X_p(ind));
            end
            Yhat = sum(F .* X_all,1);
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
                ind = ( (i-1)*obj.k_pi_per+1 : i*obj.k_pi_per );
                F(i,:) = obj.fcn_F(pi_in(ind), X_p(ind));
            end
            Yhat = sum(F .* X_all,1);
        end
        %
        function dYhat = fcn_dYhat(obj, theta_in, X_in, model_type)
            if model_type == 0
                dYhat = obj.fcn_dYhat_full(theta_in, X_in);
            elseif model_type == 1
                dYhat = obj.fcn_dYhat_pars(theta_in, X_in);
            end
        end
        function dYhat = fcn_dYhat_pars(obj, theta_in, X_in)
            X_d = X_in(obj.ind_k_delta, :);
            X_l = X_in(obj.ind_k_lambda_pars, :);
            X_p = X_in(obj.ind_kX_pi_pars, :);
            delta_in = theta_in(obj.ind_k_delta);
            lambda_in = theta_in(obj.ind_k_lambda_pars);
            pi_in = theta_in(obj.ind_k_pi_pars);
            len_l = 1;
            % g
            for i = 1:obj.k_delta
                X_d(i,:) = obj.fcn_g(X_d(i,:), 0, delta_in(i));
            end
            for i = 1:len_l
                X_l(i,:) = obj.fcn_g(X_l(i,:), 0, lambda_in(i));
            end
            g = [X_d; X_l];
            % dg
            X_d2 = X_in(obj.ind_k_delta, :);
            X_l2 = X_in(obj.ind_k_lambda_pars, :);
            for i = 1:obj.k_delta
                X_d2(i,:) = obj.fcn_dg(X_d2(i,:), 0, delta_in(i));
            end
            for i = 1:len_l
                X_l2(i,:) = obj.fcn_dg(X_l2(i,:), 0, lambda_in(i));
            end
            dg = [X_d2; X_l2];
            % F
            F = zeros(obj.k_delta+len_l, obj.n);
            for i = 1:obj.k_delta+len_l
                ind = ( (i-1)*obj.k_pi_per+1 : i*obj.k_pi_per );
                F(i,:) = obj.fcn_F(pi_in(ind), X_p(ind));
            end
            % dF
            dF = X_p;
            gdF = zeros(obj.k_pi_per*(obj.k_delta+len_l),obj.n);
            for i = 1:obj.k_delta+len_l
                ind = ( (i-1)*obj.k_pi_per+1 : i*obj.k_pi_per );
                gdF(ind,:) = dF(ind) .* repmat(g(i,:), obj.k_pi_per, 1);
            end
            %
            dYhat = [F .* dg; gdF];
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
            X_d2 = X_in(obj.ind_k_delta, :);
            X_l2 = X_in(obj.ind_k_lambda_full, :);
            for i = 1:obj.k_delta
                X_d2(i,:) = obj.fcn_dg(X_d2(i,:), 0, delta_in(i));
            end
            for i = 1:len_l
                X_l2(i,:) = obj.fcn_dg(X_l2(i,:), 0, lambda_in(i));
            end
            dg = [X_d2; X_l2];
            % F
            F = zeros(obj.k_delta+len_l, obj.n);
            for i = 1:obj.k_delta+len_l
                ind = ( (i-1)*obj.k_pi_per+1 : i*obj.k_pi_per );
                F(i,:) = obj.fcn_F(pi_in(ind), X_p(ind));
            end
            % dF
            dF = X_p;
            %
            gdF = zeros(obj.k_pi_per*(obj.k_delta+len_l),obj.n);
            for i = 1:obj.k_delta+len_l
                ind = ( (i-1)*obj.k_pi_per+1 : i*obj.k_pi_per );
                gdF(ind,:) = dF(ind) .* repmat(g(i,:), obj.k_pi_per, 1);
            end
            %
            dYhat = [F .* dg; gdF];
        end
        %
        function obj = fcn_bs2(obj, M, seed)
            % estimate null imposed model first
            Y_in = obj.Y0;
            X_null = obj.x_delta;
            bounds_pi = [-1 1]; 
            theta_hat_null_imposed = obj.fcn_estimation(0, X_null, Y_in, bounds_pi); % only the full model needs to be estimated here as this is the null imposed model
                X_d = X_null(1:obj.k_delta,:);
                delta_hat_null = theta_hat_null_imposed(1:obj.k_delta); 
                pi_hat_null = theta_hat_null_imposed(obj.k_delta+1:end);
                X_d(end,:) = obj.fcn_g(X_d(end,:), pi_hat_null(1));  
                e_hat_null = Y_in - delta_hat_null' * X_d;
%                 e_hat_null = obj.fcn_nu_hat(theta_hat_null_imposed, X_null, Y_in, 0);
            % Draw random normals
            obj.cluster_bs2_distr_seed = seed;
            rng(seed);
            Z_mat = randn(M, obj.n);
            %
            bounds_pi = [obj.LB_pi obj.UB_pi]; 
            X_all = [obj.x_delta; obj.x_lambda];
            % storage mats
            obj.theta_bs2_hat_full = zeros( obj.k_delta+1+2*obj.k_lambda, M );
            obj.theta_bs2_hat_pars = zeros( (obj.k_delta+1+1+1)*obj.k_lambda, M );
            obj.lambda_bs2_hat_full = zeros( obj.k_lambda, M );
            obj.lambda_bs2_hat_pars = zeros( obj.k_lambda, M );
            obj.max_test_stat_bs2 = zeros( M, 1 );
            obj.max_t_test_stat_bs2 = zeros( M, 1 );
            obj.wald_Cheng_test_stat_bs2 = zeros( M, 1 );
            for m = 1:M
                % Form the m'th bs model
                Z_m = Z_mat(m,:);
                temp_e = (Z_m .* e_hat_null);
                Y_m = delta_hat_null' * X_d + (Z_m .* e_hat_null);
                % Estimate the full and pars models
                temp_theta_bs_hat_full  = obj.fcn_estimation(0, X_all, Y_m, bounds_pi);
                temp_theta_bs_hat_pars  = obj.fcn_estimation(1, X_all, Y_m, bounds_pi);
                temp_lambda_bs_hat_full = temp_theta_bs_hat_full(obj.k_delta+1:obj.k_delta+obj.k_lambda);
                temp_lambda_bs_hat_pars = temp_theta_bs_hat_pars(obj.k_delta+1,1:obj.k_lambda)';
                % Store the estimates
                obj.theta_bs2_hat_full(:,m)  = temp_theta_bs_hat_full;
                obj.lambda_bs2_hat_full(:,m) = temp_lambda_bs_hat_full;
                obj.theta_bs2_hat_pars(:,m)  = temp_theta_bs_hat_pars(:);
                obj.lambda_bs2_hat_pars(:,m) = temp_lambda_bs_hat_pars;
                % form the test stats
                theta_in = temp_theta_bs_hat_full;
                [~, V_hat_full_lambda] = obj.fcn_V_hat(theta_in, X_all, Y_m, 0, 0);
                theta_in = temp_theta_bs_hat_pars;
                [~, V_hat_pars_lambda] = obj.fcn_V_hat(theta_in, X_all, Y_m, 1, 0);
                % Cheng's Wald Stat
                obj.wald_Cheng_test_stat_bs2(m) = obj.n * (temp_lambda_bs_hat_full)' * inv(V_hat_full_lambda) * (temp_lambda_bs_hat_full);
                % Max Test Stats
                W = ones(obj.k_lambda,1);
                Wt = sqrt(diag(V_hat_pars_lambda)).^(-1);
                obj.max_test_stat_bs2(m) = sqrt(obj.n) * max(abs(W .* (temp_lambda_bs_hat_pars)));
                obj.max_t_test_stat_bs2(m) = sqrt(obj.n) * max(abs(Wt .* (temp_lambda_bs_hat_pars)));
            end
        end
        %
        function out = fcn_F(obj, pi, X_pi)
            out = pi' * X_pi;
        end
        function out = fcn_g(obj, X2, c, gamma) %g_E2 swaps c and pi %fcn_g_E2
            out =  1 - exp( -gamma .* (X2 - c).^2 );
        end
        function out = fcn_dg(obj, X2, c, gamma)
            out = ( ((X2 - c).^2) .* exp( -gamma .* (X2 - c).^2 ) );
        end
        function out = fcn_d2g(obj, X2, c, gamma)
            out = -( ((X2 - c).^4) .* exp( -gamma .* (X2 - c).^2 ) );
        end
        %
        function nu_hat = fcn_nu_hat(obj, theta, X, Y, type, Taylor_flag)
            if nargin < 6
                Taylor_flag = 0;
            end
            if type == 0
                nu_hat = obj.fcn_nu_hat_full(theta, X, Y, Taylor_flag);
            elseif type == 1
                nu_hat = obj.fcn_nu_hat_pars(theta, X, Y, Taylor_flag);
            end
        end
        function nu_hat = fcn_nu_hat_full(obj, theta, X, Y, Taylor_flag)
            Yhat = obj.fcn_Yhat(theta, X, 0);                %
            nu_hat = Y - Yhat;
        end
        function nu_hat = fcn_nu_hat_pars(obj, theta, X, Y, Taylor_flag)
            Yhat = obj.fcn_Yhat(theta, X, 1);                %
            nu_hat = Y - Yhat;
        end
        %
        function D_nu_hat = fcn_D_nu_hat(obj, theta, X, type, type_exclude_beta, Taylor_flag)
            if nargin < 5
                type_exclude_beta = 0;
            end
            if nargin < 6
                Taylor_flag = 0;
            end
            if type == 0
                D_nu_hat = obj.fcn_D_nu_hat_full(theta, X, type_exclude_beta);
            elseif type == 1
                D_nu_hat = obj.fcn_D_nu_hat_pars(theta, X, type_exclude_beta);
            end
        end        
        function D_nu_hat = fcn_D_nu_hat_full(obj, theta, X, type_exclude_beta, Taylor_flag)
            D_nu_hat = -obj.fcn_dYhat(theta, X, 0);
        end
        function D_nu_hat = fcn_D_nu_hat_pars(obj, theta, X, type_exclude_beta, Taylor_flag)
            D_nu_hat = -obj.fcn_dYhat(theta, X, 1);
        end
        %
        function m = fcn_m(obj, theta, X, Y, type)
            nu_hat = obj.fcn_nu_hat(theta, X, Y, type);
            m = nu_hat.^2 / 2;
        end
        function Dm = fcn_Dm(obj, theta, X, Y, type, type_exclude_beta, nu_hat)
            if nargin < 6
                type_exclude_beta = 0;
            end
            if nargin < 7
                nu_hat = obj.fcn_nu_hat(theta, X, Y, type);
            end
            D_nu_hat = obj.fcn_D_nu_hat(theta, X, type, type_exclude_beta);
            Dm = repmat(nu_hat, size(D_nu_hat,1),1) .* D_nu_hat;
        end
        %
        function l_K_beta = fcn_l_K_beta(obj, type, l_K)
            if type == 0 % full model so all lambdas
                % l_K_beta = l_K - (obj.k_delta + obj.k_lambda);
                ind_beta = [obj.ind_k_delta obj.ind_k_lambda_full];
                ind_pi = obj.ind_k_pi_full;
                % [sharedvals, idx] = intersect(A,B,'stable');
                % [tf,idx] = ismember(B,A)
            elseif type == 1 % Pars model so only 1 lambda
                % l_K_beta = l_K - (obj.k_delta + 1);
                ind_beta = [obj.ind_k_delta obj.ind_k_lambda_pars];
                ind_pi = obj.ind_k_pi_pars;
            end
            [~, i2] = intersect(ind_pi, l_K, 'stable');
            C = kron(ind_beta, ones(1,obj.k_pi_per));
            l_K_beta = unique(C(i2));
        end
        function G = fcn_G(obj, theta, X, Y, type, l_K, type_exclude_beta, Z)
            % l_K_beta = obj.fcn_l_K_beta(type, l_K);
            % theta(l_K_beta) = zeros(length(l_K_beta),1);
            method = 2;
            if method == 1
                if nargin < 8
                    Dm = obj.fcn_Dm(theta, X, Y, type, type_exclude_beta);
                    Dm = Dm - (sum(Dm,2) / obj.n);
                elseif nargin == 8
                    Dm = obj.fcn_Dm(theta, X, Y, type, type_exclude_beta, Z);
                end
                Dm(l_K,:) = [];
            elseif method == 2
                if nargin < 8
                    Z = ones(1, obj.n);
                end
                Dm = obj.fcn_Dm(theta, X, Y, type, type_exclude_beta);
                Dm = Dm - (sum(Dm,2) / obj.n);
                Dm(l_K,:) = [];
                Dm = Dm .* repmat(Z, size(Dm,1), 1);
            end
            G = Dm / sqrt(obj.n);
        end
        function H = fcn_H(obj, theta, X, Y, type, l_K, type_exclude_beta)
            % l_K_beta = obj.fcn_l_K_beta(type, l_K);
            % theta(l_K_beta) = zeros(length(l_K_beta),1);
            D_nu_hat = obj.fcn_D_nu_hat(theta, X, type, type_exclude_beta);
            D_nu_hat(l_K,:) = [];
            H = D_nu_hat * D_nu_hat' / obj.n;
        end
        function K = fcn_K(obj, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in)
            % l_K_beta = obj.fcn_l_K_beta(type, l_K);
            % theta(l_K_beta) = zeros(length(l_K_beta),1);
            D_nu_hat = obj.fcn_D_nu_hat(theta, X, type, type_exclude_beta);
            theta2 = theta; theta2(l_K) = theta_0_in(l_K);  % recall l_K is the index location of pi_{l_K}, not beta
            D_nu_hat_0 = obj.fcn_D_nu_hat(theta2, X, type, type_exclude_beta);
            D_nu_hat(l_K,:) = [];
            % % l_K indexes the weakly id'd pi.  We must locate the
            % % associated beta's.  The size of theta is 
            % % k_delta + k_lambda + 1 + k_lambda.  The betas are
            % % [k_delta(end) k_delta+1:k_delta+k_lambda].
            % l_K_beta = l_K - (obj.k_lambda + 1);
            % S = zeros(1, length(theta)); S(l_K_beta) = ones(1, length(l_K_beta));
            % temp = S * D_nu_hat_0;
            % I do this outside of the function now.
            temp = D_nu_hat_0;
            K = D_nu_hat * temp' / obj.n;
        end
        %
        function tau = fcn_tau(obj, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z)
            l_K_beta = obj.fcn_l_K_beta(type, l_K);
            theta(l_K_beta) = zeros(length(l_K_beta),1);
            G = obj.fcn_G(theta, X, Y, type, l_K, type_exclude_beta, Z);
            G = sum(G,2);
            H = obj.fcn_H(theta, X, Y, type, l_K, type_exclude_beta);
            K = obj.fcn_K(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in);
            % l_K_beta = obj.fcn_l_K_beta(type, l_K);
            % S = zeros(1, length(theta)); S(l_K_beta) = ones(1, length(l_K_beta));
            % b_l_K(l_K_beta) = theta_0_in(l_K_beta) * sqrt(obj.n);
            Sb = zeros(1, length(theta))';
            Sb(l_K_beta) = theta_0_in(l_K_beta) * sqrt(obj.n);
            if isempty(l_K) ~= 1
                % tau = inv(H) * ( K * b_l_K + G );
                % tau = inv(H) * ( K * Sb + G );
                tau = inv(H) * ( G - K * Sb );
            elseif isempty(l_K) == 1
                tau = inv(H) * G;
            end
        end
        function chi = fcn_chi(obj, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z)
            l_K_beta = obj.fcn_l_K_beta(type, l_K);
            theta(l_K_beta) = zeros(length(l_K_beta),1);
            H = obj.fcn_H(theta, X, Y, type, l_K, type_exclude_beta);
            tau = obj.fcn_tau(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
            chi = - (1/2) * tau' * H * tau;
        end
        function chi = fcn_pi_star_objective(obj, pi_in, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z)
            theta(l_K) = pi_in;
            chi = obj.fcn_chi(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
        end
        function pi_star = fcn_pi_star(obj, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z)
            ind = obj.fcn_l_K_beta(type, l_K);
            LB_pi_K = obj.LB_pi(ind);
            UB_pi_K = obj.UB_pi(ind);
            step_size = .1;
            % chi_temp = inf;
            % First get a starting point
            pi_init = (UB_pi_K - LB_pi_K) .* rand(size(LB_pi_K,1),size(LB_pi_K,2)) + LB_pi_K;
            chi_temp = obj.fcn_pi_star_objective(pi_init, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
            for iter = 1:5 % iterate the line search a few times
            for ind2 = 1:length(ind) % simple line search on the objective function
                pi_temp = pi_init;
                for temp = (LB_pi_K(ind2):step_size:UB_pi_K(ind2))
                    pi_temp(ind2) = temp;
                    chi_temp2 = obj.fcn_pi_star_objective(pi_temp, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
                    % fprintf('%e %e %.4f %.4f \n', chi_temp2, chi_temp, pi_init(ind2), pi_temp(ind2))
                    if (chi_temp2 < chi_temp) == 1
                        chi_temp = chi_temp2;
                        pi_init(ind2) = pi_temp(ind2);
%                         fprintf('%.4f \n', pi_init);
                    end
                end
            end
            end
            % Now use a gradient descent method:
            [pi_star, fval_temp, exit_flag_temp] = fmincon(@(pi_in) obj.fcn_pi_star_objective(pi_in, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z), pi_init, [], [], [], [], LB_pi_K, UB_pi_K, [], obj.options_fmincon);
            if exit_flag_temp < 0  % for error checking
                pi_star = NaN;
            end
            if (chi_temp < fval_temp) == 1 % Not sure why the gradient descent misses this sometimes...
                pi_star = pi_init;
            end            
        end
        function [mathfrac_z, HG] = fcn_mathfrac_z(obj, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z)
            l_K_beta = obj.fcn_l_K_beta(type, l_K);
            theta_in_temp = theta;
            theta(l_K_beta) = zeros(length(l_K_beta),1);
            pi_star = obj.fcn_pi_star(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
            theta(l_K) = pi_star;
            tau = fcn_tau(obj, theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
            % l_K_beta = obj.fcn_l_K_beta(type, l_K);
            b_l_K = zeros(size(tau,1),size(tau,2));
            b_l_K(l_K_beta) = theta_0_in(l_K_beta) * sqrt(obj.n);
            % tau(l_K_beta) = tau(l_K_beta) - b_l_K;
            tau = tau - b_l_K;
            mathfrac_z = zeros(length(theta), 1);
            ind_all = (1:length(theta));
            ind_not_l_K = setdiff(ind_all, l_K);
            mathfrac_z(ind_not_l_K) = tau;
%             mathfrac_z(l_K) = abs(mathfrac_z(l_K_beta)) .* pi_star;  % oops!
            mathfrac_z(l_K) = abs(mathfrac_z(l_K_beta) + b_l_K(l_K_beta)) .* pi_star;
%             theta = obj.theta_0;
            theta = theta_in_temp;
            theta(l_K) = pi_star;
            G = obj.fcn_G(theta, X, Y, type, [], type_exclude_beta, Z);  
            H = obj.fcn_H(theta, X, Y, type, [], type_exclude_beta); 
            HG = inv(H) * G;
        end
        function [mathfrac_z, V]= fcn_mathfrac_z_wrapper(obj, type, Z)
            Y = obj.Y0;
            type_exclude_beta = 1;  %%%%%%%%%%%%%%%%              
            theta_temp = obj.theta_0;
            theta_temp(obj.ind_k_lambda_full) = zeros(obj.k_lambda,1);
            X_in = [obj.x_delta; obj.x_lambda; obj.x_pi];
            len_full = length(obj.theta_hat_full);
            len_pars = length(obj.theta_hat_pars(:));
            if type == 0
                X = X_in;
                theta = obj.theta_hat_full;        
%                 theta = obj.theta_0;
%                 theta_0_in = obj.theta_0;
%                 theta = theta_temp;
                theta_0_in = theta_temp;
                i = [];
                l_K = obj.fcn_l_K_selection(theta, X, Y, type, i, type_exclude_beta);
                [mathfrac_z, HG]= obj.fcn_mathfrac_z(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
                V = HG * HG';
            elseif type == 1
                mathfrac_z = zeros(len_pars,1);
                HG = zeros(len_pars, obj.n);
                for i = 1:obj.k_lambda
                    ind_pi = [(1:obj.k_pi_per) ( (i-1+1)*obj.k_pi_per+1 : (i+1)*obj.k_pi_per )];
                    ind_theta_i = [(1:obj.k_delta) obj.k_delta+i obj.k_delta+obj.k_lambda+ind_pi];
                    X = X_in(ind_theta_i, :);
                    theta = obj.theta_hat_pars(:,i);        
%                     theta = obj.theta_0(ind_theta_i);
%                     theta_0_in = obj.theta_0(ind_theta_i);
%                     theta = theta_temp(ind_theta_i);
                    theta_0_in = theta_temp(ind_theta_i);
                    l_K = obj.fcn_l_K_selection(theta, X, Y, type, i, type_exclude_beta);
                    temp = (obj.k_delta+1+obj.k_pi_per+obj.k_pi_per);
                    ind = ( temp*(i-1)+1:temp*(i) );
                    [mathfrac_z(ind), HG(ind,:)]= obj.fcn_mathfrac_z(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
                end
                V = HG * HG';
            end                
        end
        function [mathfrac_z_distr, V] = fcn_mathfrac_z_distr(obj, type, M, Z_mat)
            % M = 100;
            if type == 0
                len = length(obj.theta_hat_full);
                mathfrac_z_distr = zeros(len, M);
                V = zeros(len, len, M);
            elseif type == 1
                len = length(obj.theta_hat_pars(:));
                mathfrac_z_distr = zeros(len, M);
                V = zeros(len, len, M);
            end
            % Z_mat = randn(M, obj.n);
            for m = 1:M
                Z = Z_mat(m,:);
                [temp, temp_V] = obj.fcn_mathfrac_z_wrapper(type, Z);
                mathfrac_z_distr(:,m) = temp;
                V(:,:,m) = temp_V;
            end
        end
        function obj = fcn_cluster_sim_distr(obj, M, seed)
%             M = 10; % simulate on cluster, so these will be combined ex post.
%             seed = obj.seed_input;
            obj.cluster_sim_distr_seed = seed;
            rng(seed);
            Z_mat = randn(M, obj.n);
            type = 0;
            [obj.mathfrac_z_distr_full, V_full]= obj.fcn_mathfrac_z_distr(type, M, Z_mat);
            type = 1;
            [obj.mathfrac_z_distr_pars, V_pars] = obj.fcn_mathfrac_z_distr(type, M, Z_mat);
            lambda_ind_full = (obj.k_delta+1:obj.k_delta+obj.k_lambda);
            obj.lambda_distr_full = obj.mathfrac_z_distr_full(lambda_ind_full,:);
            lambda_ind_pars = ( obj.k_delta+1 : (obj.k_delta+1+1+1) : (obj.k_delta+1+2)*(obj.k_lambda) );
            obj.lambda_distr_pars = obj.mathfrac_z_distr_pars(lambda_ind_pars,:);
            V_lam_full = V_full(lambda_ind_full, lambda_ind_full, :);
            V_lam_pars = V_pars(lambda_ind_pars, lambda_ind_pars, :);
            obj.test_distr_wald_Cheng = cell(obj.k_lambda,1);
            obj.test_distr_max = cell(obj.k_lambda,1);
            obj.test_distr_max_t = cell(obj.k_lambda,1);
            for k = 1:obj.k_lambda
                obj.test_distr_wald_Cheng{k} = zeros(M,1);
                obj.test_distr_max{k} = zeros(M,1);
                obj.test_distr_max_t{k} = zeros(M,1);
                for m = 1:M
                    lambda_full_m = obj.lambda_distr_full(:,m);
                    V_lam_full_m = V_lam_full(:,:,m);
                    lambda_pars_m = obj.lambda_distr_pars(:,m);
                    V_lam_pars_m = V_lam_pars(:,:,m);
                    Wt = 1 ./ sqrt(diag(V_lam_pars_m));
                    obj.test_distr_wald_Cheng{k}(m) = lambda_full_m(1:k)' * inv(V_lam_full_m(1:k,1:k)) * lambda_full_m(1:k);
                    obj.test_distr_max{k}(m) = max(abs(lambda_pars_m(1:k)));
                    obj.test_distr_max_t{k}(m) = max(abs(Wt(1:k) .* lambda_pars_m(1:k)));
                end
            end
        end
        %
        function l_K_out = fcn_l_K_selection(obj, theta, X, Y, type, i, type_exclude_beta)
            % first get the beta indices:
            if type == 0 % full
                % ind_beta = (obj.k_delta:obj.k_delta+obj.k_lambda);
                % ind_pi = ind_beta + (obj.k_lambda + 1);
                ind_beta = [obj.ind_k_delta obj.ind_k_lambda_full];
                ind_pi = obj.ind_k_pi_full;
                theta_temp = obj.theta_hat_full;
            elseif type == 1 % pars
                % ind_beta = [obj.k_delta obj.k_delta+1];
                % ind_pi = ind_beta + (1 + 1);
                ind_beta = [obj.ind_k_delta obj.ind_k_lambda_pars];
                ind_pi = obj.ind_k_pi_pars;
                theta_temp = obj.theta_hat_pars(:,i);
            end
            l_K_in = [];
            Z = ones(1, obj.n);
            G = obj.fcn_G(theta_temp, X, Y, type, l_K_in, type_exclude_beta, Z);
            H = obj.fcn_H(theta_temp, X, Y, type, l_K_in, type_exclude_beta);
%             Sigma = inv(H) * (G * G') * inv(H)' / sqrt(obj.n);
            Sigma = inv(H) * (G * G') * inv(H)';
%             inv_Sigma = inv(Sigma);
            l_K_beta = [];
            for j = 1 % for this experiment, H_0: beta = lambda = 0, so we put all the lambdas in l_K
                % in other experiments, we use this procedure on all betas
                beta_j = theta_temp(ind_beta(j));
                ICS = sqrt( obj.n * beta_j * ( 1/Sigma(j,j) ) * beta_j );
                if (ICS <= obj.kappa_n) == 1
                    l_K_beta = [l_K_beta ind_beta(j)]; % include j's associated pi
                elseif (ICS > obj.kappa_n) == 1
                    % exclude j: do nothing here
                end
            end
            l_K_beta = [l_K_beta ind_beta(2:end)]; % for this experiment, H_0: beta = lambda = 0, so we put all the lambdas in l_K
            C = kron(ind_beta, ones(1,obj.k_pi_per));
            l_K_out = ind_pi(find(ismember(C, l_K_beta)));
        end
        %
        function obj = fcn_estimation_initial(obj)
            Y_in = obj.Y0;
            X_all = [obj.x_delta; obj.x_lambda; obj.x_pi];
            obj.LB_pi = -2*ones(obj.k_pi, 1);
            obj.UB_pi = -obj.LB_pi;
            bounds_pi = [obj.LB_pi obj.UB_pi];
            bnd_dl = repmat(3,obj.k_delta+obj.k_lambda,1);
            LB = [-bnd_dl; bounds_pi(:,1)]; 
            UB = [ bnd_dl; bounds_pi(:,2)];
            bounds_in = [LB UB];
            theta_init = [obj.delta_0; obj.lambda_0; obj.pi_0];
            c = 2;
            theta_init = theta_init -c/2 + c*rand(length(theta_init),1);
            obj.theta_hat_full = obj.fcn_estimation(0, X_all, Y_in, bounds_in, theta_init);
            obj.theta_hat_pars = obj.fcn_estimation(1, X_all, Y_in, bounds_in, theta_init);
            %
            obj.lambda_hat_full = obj.theta_hat_full(obj.ind_k_lambda_full);
            obj.lambda_hat_pars = obj.theta_hat_pars(obj.ind_k_lambda_pars,1:obj.k_lambda)';
            % [lambda_hat_full lambda_hat_pars]
        end
        % Estimation Function below
        function theta_out = fcn_estimation(obj, type, X_in, Y_in, bounds, theta_init)
            if type == 0
                theta_out = obj.fcn_estimation_full_2(type, X_in, Y_in, bounds, theta_init);
            elseif type == 1
                theta_out = obj.fcn_estimation_pars_2(type, X_in, Y_in, bounds, theta_init);
            end
        end
        % Estimation via fmincon here
        function [Q, DQ] = fcn_loss(obj, theta, X, Y, type)
            m = obj.fcn_m(theta, X, Y, type);
            Dm = obj.fcn_Dm(theta, X, Y, type);
            Q = sum(m) / obj.n;
            DQ = sum(Dm,2) / obj.n;
            DQ = DQ';
        end
        function theta_out = fcn_estimation_full_2(obj, type, X_in, Y_in, bounds_in, theta_init)
            LB = bounds_in(:,1);
            UB = bounds_in(:,2);
            %             bnd_p = 1/sqrt(obj.n) .* [1; zeros(size(bounds_pi,1)-1,1)];
            %             bnd_dl = repmat(3,obj.k_delta+obj.k_lambda,1);
            %             LB = [-bnd_dl; bounds_pi(:,1)]; 
            %             UB = [bnd_dl; -bnd_p];
            %             [theta_out, fval_temp, exit_flag_temp] = fmincon(@(theta_in) obj.fcn_loss(theta_in, X_in, Y_in, type), theta_init, [], [], [], [], LB, UB, [], obj.options_fmincon);
            %             LB = [-bnd_dl; bnd_p]; 
            %             UB = [bnd_dl; bounds_pi(:,2)];
            %             rng(obj.rng_seed);
            %             theta_init = LB + (UB-LB).* rand(length(LB),1);
%             [theta_out, fval_temp2, exit_flag_temp] = fmincon(@(theta_in) obj.fcn_loss(theta_in, X_in, Y_in, type), theta_init, [], [], [], [], LB, UB, [], obj.options_fmincon);
            [theta_out, fval_temp2, exit_flag_temp] = lsqnonlin(@(theta_in) obj.fcn_loss(theta_in, X_in, Y_in, type), theta_init, [], [], obj.options_lsqnonlin);
            % lsqnonlin(fun,x0,lb,ub,options)
            %             if fval_temp2 < fval_temp
            %                 theta_out = theta_out2;
            %             end
        end
        function theta_out = fcn_estimation_pars_2(obj, type, X_in, Y_in, bounds_in, theta_init_in)
            theta_out = zeros((obj.k_delta+1)*(obj.k_pi_per+1), obj.k_lambda);
            LB_all = bounds_in(:,1);
            UB_all = bounds_in(:,2);
            %             bnd_p = 1/sqrt(obj.n) .* [1; zeros(size(bounds_pi,1)-1,1)];
            %             bnd_dl = repmat(3,obj.k_delta+obj.k_lambda,1);
            %             LB_all1 = [-bnd_dl; bounds_pi(:,1)]; 
            %             UB_all1 = [bnd_dl; -bnd_p];
            %             LB_all2 = [-bnd_dl; bnd_p]; 
            %             UB_all2 = [bnd_dl; bounds_pi(:,2)];
            %             X_d = X_in(obj.ind_k_delta, :);
            %             X_l = X_in(obj.ind_k_lambda_full, :);
            %             X_p = X_in(obj.ind_kX_pi_full, :);
            for i = 1:obj.k_lambda
                ind_pi = [(1:obj.k_pi_per) ((i*obj.k_pi_per+1):((i+1)*obj.k_pi_per))]; 
                ind_pars = [1:obj.k_delta obj.k_delta+i obj.k_delta+obj.k_lambda+ind_pi];
                theta_init = theta_init_in(ind_pars);
                X_pars = X_in(ind_pars,:);
                LB = LB_all(ind_pars); 
                UB = UB_all(ind_pars);
                %                 LB = LB_all1(ind_pars); 
                %                 UB = UB_all1(ind_pars);
                %                 [theta_out_temp, fval_temp, exit_flag_temp] = fmincon(@(theta_in) obj.fcn_loss(theta_in, X_pars, Y_in, type), theta_init, [], [], [], [], LB, UB, [], obj.options_fmincon);
                %                 LB = LB_all2(ind_pars); 
                %                 UB = UB_all2(ind_pars);
                %                 rng(obj.rng_seed);
                %                 theta_init = LB + (UB-LB).* rand(length(LB),1);
%                 [theta_out_temp, fval_temp2, exit_flag_temp] = fmincon(@(theta_in) obj.fcn_loss(theta_in, X_pars, Y_in, type), theta_init, [], [], [], [], LB, UB, [], obj.options_fmincon);
                [theta_out_temp, fval_temp2, exit_flag_temp] = lsqnonlin(@(theta_in) obj.fcn_loss(theta_in, X_pars, Y_in, type), theta_init, [], [], obj.options_lsqnonlin);
                %                 if fval_temp2 < fval_temp
                %                     theta_out_temp = theta_out_temp2;
                %                 end
                theta_out(:,i) = theta_out_temp;
            end
        end
        %
        function [B_full, B_pars] = fcn_B(obj, theta_full, theta_pars)
            B_full = [ones(obj.k_delta+obj.k_lambda,1); theta_full(obj.k_delta:obj.k_delta+obj.k_lambda)];
            num = (obj.k_delta+1+1+1);
            B_pars = ones( num*obj.k_lambda,1);
            ind = [];
            for temp_k = 1:obj.k_lambda
                temp = num*(temp_k-1) + (obj.k_delta+1+1:num);
                ind = [ind temp];
            end
            ind2 = ind-(obj.k_delta+1);
            B_pars(ind) = theta_pars(ind2);
        end
        %
    end % end methods
end % end class

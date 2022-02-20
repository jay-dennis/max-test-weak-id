classdef class_dgp_2
    properties (SetAccess=public, GetAccess=public)
	 % dgp variables
     test_number
     dgp_type, dgp_type_string
     hypothesis_type, hypothesis_type_string
     g_type, g_type_string
     n, init_n
	 Y0
     e
     k_delta, k_lambda
     x_delta, x_beta2, x_pi
     zeta_0, beta_0, pi_0, theta_0
     delta_0, beta2_0
     b1_0, b2_0
     ind_lambda_full, ind_lambda_pars
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
     options_fmincon     
    end
    methods
        function obj = class_dgp_1(dgp_type, hypothesis_type, n, seed_input, k_lambda, k_delta, beta_in)
            obj.test_number = 1;
            % Define the location of the lambdas:
            obj.ind_lambda_full = (obj.k_delta+1:obj.k_delta+obj.k_lambda);  % beta2 location
            obj.ind_lambda_pars = obj.k_delta+1;  % beta2 location
            % obj.ind_lambda_pars_all_stacked = ( obj.k_delta+1 : (obj.k_delta+1+1+1) : (obj.k_delta+1+2)*(obj.k_lambda) ); % not needed I think
            %
            obj.g_type = 2; 
                if obj.g_type == 1 
                    obj.g_type_string = 'Logistic';
                elseif obj.g_type == 2
                    obj.g_type_string = 'Exponential';
                end
            obj.dgp_type = dgp_type;
            obj.hypothesis_type = hypothesis_type;
            obj.n = n;
            obj.k_delta = k_delta; 
            obj.k_lambda = k_lambda;
            obj.seed_input = seed_input;
            obj.kappa_n = sqrt(log(obj.n)); % for ICS cvs
            obj.options_fmincon = optimoptions(@fmincon,'Display','off','Diagnostics','off'); 
            %
            % Now generate the data
	        obj.init_n = obj.n + 0;
            rng(0); % reproduce the X's
            obj.x_pi = randn(1,obj.n); % the transition variable
            switch dgp_type
                case 1
                    obj.dgp_type_string = 'Case 1: iid';
                    obj.x_delta = randn(obj.k_delta,obj.n);
                    obj.x_beta2 = randn(obj.k_lambda,obj.n);
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
                    obj.x_beta2 = Sig * w_beta;
                case 3 
                    obj.dgp_type_string = 'Case 3: with-in and cross-block dependence';
                    num = obj.k_delta + obj.k_lambda;
                    temp = randn(num,obj.n);
                    rho = .5;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num));
                    Sig = chol(Sig);
                    X = Sig * temp;
                    obj.x_delta = X(1:obj.k_delta,:);
                    obj.x_beta2 = X(obj.k_delta+1:end,:);
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
            
            obj.beta2_0 = beta_in(2) * ones(obj.k_lambda,1);

            obj.beta_0 = beta_in(1);
            obj.delta_0 = ones(obj.k_delta,1);
            obj.delta_0(end) = beta_in(1);
            obj.pi_0 = zeros(obj.k_lambda+1,1);

            obj.theta_0 = [obj.delta_0; obj.beta2_0; obj.pi_0];
            
            obj.b1_0 = beta_in(1) * sqrt(obj.n);
            obj.b2_0 = obj.beta2_0 * sqrt(obj.n);
            
            if obj.k_delta > 1 % Note: Make sure that the inputed k_delta > 1
                obj.delta_0(1) = 0; % intercept = 0 in true model.
                obj.x_delta(1,:) = ones(1, obj.n); % replace the 1st variable with the intercept.
            end
            
            X_d = obj.x_delta;
            X_d(end,:) = obj.fcn_g(obj.x_delta(end,:), obj.pi_0(1));
            
            X_b2 = obj.x_beta2;
            for i = 1:obj.k_lambda
                X_b2(i,:) = obj.fcn_g(obj.x_beta2(i,:), obj.pi_0(i+1));
            end
            
            obj.Y0 = obj.delta_0' * X_d + obj.beta2_0' * X_b2 + obj.e;
            
            obj.Y0 = obj.Y0((obj.init_n - obj.n + 1):end); % remove Y0 burn-in values
        end
        function obj = fcn_estimation_initial(obj)
            Y_in = obj.Y0;
            X_all = [obj.x_delta; obj.x_beta2];
%             % X_pi = [obj.x_delta(end,:); obj.x_beta2];
%             X_pi = repmat(obj.x_pi, 1+obj.k_lambda, 1);
%             X_pi_sorted = sort(X_pi,2);
%             obj.LB_pi = X_pi_sorted(:,floor(.15 * obj.n));
%             obj.UB_pi = X_pi_sorted(:,floor(.85 * obj.n)); 
            obj.LB_pi = -1*ones(1+obj.k_lambda, 1);
            obj.UB_pi = ones(1+obj.k_lambda, 1);
            bounds_pi = [obj.LB_pi obj.UB_pi]; 
            obj.theta_hat_full = obj.fcn_estimation(0, X_all, Y_in, bounds_pi);
            obj.theta_hat_pars = obj.fcn_estimation(1, X_all, Y_in, bounds_pi);
            %
            obj.lambda_hat_full = obj.theta_hat_full(obj.ind_lambda_full);
            obj.lambda_hat_pars = obj.theta_hat_pars(obj.ind_lambda_pars,1:obj.k_lambda)';
            % [lambda_hat_full lambda_hat_pars]
        end
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
            X_all = [obj.x_delta; obj.x_beta2];
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
                temp_lambda_bs_hat_full = temp_theta_bs_hat_full(obj.ind_lambda_full);
                temp_lambda_bs_hat_pars = temp_theta_bs_hat_pars(obj.ind_lambda_pars,1:obj.k_lambda)';
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
        function out = fcn_g(obj, X, pi)
            c = 10; X1 = X; X2 = obj.x_pi;
            if obj.g_type == 1
                out = obj.fcn_g_L(X1, X2, c, pi);
            elseif obj.g_type == 2
                out = obj.fcn_g_E(X1, X2, c, pi);
            elseif obj.g_type == 3
                out = obj.fcn_g_E2(X1, X2, c, pi);
            end
        end
        function out = fcn_dg(obj, X, pi)
            c = 10; X1 = X; X2 = obj.x_pi;
            if obj.g_type == 1
                out = obj.fcn_dg_L(X1, X2, c, pi);
            elseif obj.g_type == 2
                out = obj.fcn_dg_E(X1, X2, c, pi);
            elseif obj.g_type == 3
                out = obj.fcn_dg_E2(X1, X2, c, pi);
            end
        end
        function out = fcn_d2g(obj, X, pi)
            c = 10; X1 = X; X2 = obj.x_pi;
            if obj.g_type == 1
                out = obj.fcn_d2g_L(X1, X2, c, pi);
            elseif obj.g_type == 2
                out = obj.fcn_d2g_E(X1, X2, c, pi);
            elseif obj.g_type == 3
                out = obj.fcn_d2g_E2(X1, X2, c, pi);
            end
        end
        function out = fcn_g_L(obj, X1, X2, c, pi)
            out = X1 ./ ( 1 - exp(-c * (X2 - pi)) );
        end
        function out = fcn_dg_L(obj, X1, X2, c, pi)
            out = X1 .* ( 1 - exp(-c * (X2 - pi)) ).^(-2) .* exp(-c * (X2 - pi)) .* c;
        end
        function out = fcn_d2g_L(obj, X1, X2, c, pi)
            out = ( X1 .* 2 .* ( 1 - exp(-c * (X2 - pi)) ).^(-3) .* exp(-c * (X2 - pi)).^(2) .* (c^2) ) + ( X1 .* ( 1 - exp(-c * (X2 - pi)) )^(-2) .* exp(-c * (X2 - pi)) .* (c^2) );
        end
        %
        function out = fcn_g_E(obj, X1, X2, c, pi)
            out = X1 .* ( 1 - exp( -c .* (X2 - pi).^2 ) );
        end
        function out = fcn_dg_E(obj, X1, X2, c, pi)
            out = X1 .* ( -2 .* c .* (X2 - pi) .* exp( -c .* (X2 - pi).^2 ) );
        end
        function out = fcn_d2g_E(obj, X1, X2, c, pi)
            out = X1 .* 2 .* c .* exp( -c .* (X2 - pi).^2 ) .* ( 1 - 2 .* c .* (X2 - pi).^2 );
        end
        %
        function out = fcn_g_E2(obj, X1, X2, c, pi) %g2 swaps c and pi
            out = X1 .* ( 1 - exp( -pi .* (X2 - c).^2 ) );
        end
        function out = fcn_dg_E2(obj, X1, X2, c, pi)
            out = X1 .* ( ((X2 - pi).^2) .* exp( -c .* (X2 - pi).^2 ) );
        end
        function out = fcn_d2g_E2(obj, X1, X2, c, pi)
            out = 0; % placeholder
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
            num_var = size(X,1);
            if num_var == obj.k_delta % for the null imposed model
                X_d = X(1:obj.k_delta,:); X_b2 = 0;
                delta = theta(1:obj.k_delta); beta2 = 0;
                pi = theta(obj.k_delta+1);
                %
                X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
                %
                nu_hat = Y - delta' * X_d - beta2' * X_b2;
            elseif num_var > obj.k_delta
                X_d = X(1:obj.k_delta,:); X_b2 = X(obj.k_delta+1:end,:);
                delta = theta(1:obj.k_delta); beta2 = theta(obj.k_delta+1:obj.k_delta+obj.k_lambda);
                pi = theta(obj.k_delta+obj.k_lambda+1:end);
                %
                X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
                for i = 1:obj.k_lambda
                    X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
                end
                %
                nu_hat = Y - delta' * X_d - beta2' * X_b2;
            end
        end
%         function nu_hat = fcn_nu_hat_full(obj, theta, X, Y)
% %             X_d = X(1:obj.k_delta,:); X_b2 = X(obj.k_delta+1:end,:);
% %             delta = theta(1:obj.k_delta); beta2 = theta(obj.k_delta+1:obj.k_delta+obj.k_lambda);
% %             pi = theta(obj.k_delta+obj.k_lambda+1:end);
% %             %
% %             X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
% %             for i = 1:obj.k_lambda
% %                 X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
% %             end
%             %
%             num_var = size(X,1);
%             X_d = X(1:obj.k_delta,:);
%             delta = theta(1:obj.k_delta);
%             if num_var == obj.k_delta  % only doing this for the full model b/c the null imposed model is only estimated with the full model estimation procedure
%                 pi = theta(obj.k_delta+1:end);
%                 X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));
%                 nu_hat = Y - delta' * X_d;
%             elseif num_var > obj.k_delta
%                 X_b2 = X(obj.k_delta+1:end,:);
%                 beta2 = theta(obj.k_delta+1:obj.k_delta+obj.k_lambda);
%                 pi = theta(obj.k_delta+obj.k_lambda+1:end);
%                 X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));
%                 for i = 1:obj.k_lambda
%                     X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
%                 end
%                 nu_hat = Y - delta' * X_d - beta2' * X_b2;
%             end
%         end
        function nu_hat = fcn_nu_hat_pars(obj, theta, X, Y, Taylor_flag)
            X_d = X(1:obj.k_delta,:); X_b2 = X(obj.k_delta+1,:);
            delta = theta(1:obj.k_delta); beta2 = theta(obj.k_delta+1);
            pi = theta(obj.k_delta+2:end);
            %
            X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
            i = 1;
            X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
            %
            nu_hat = Y - delta' * X_d - beta2' * X_b2;
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
            num_var = size(X,1);
            if num_var == obj.k_delta
                X_d = X(1:obj.k_delta,:); X_b2 = 0;
                delta = theta(1:obj.k_delta); beta2 = 0;
                pi = theta(obj.k_delta+1);
                %
                X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
                %
                temp1 = obj.fcn_dg(X(obj.k_delta,:), pi(1));
                %
                if type_exclude_beta == 0
                    temp1 = delta(end) * temp1;
                elseif type_exclude_beta == 1
                    %do nothing
                end
                D_nu_hat = -[X_d; temp1];
            elseif num_var > obj.k_delta
                X_d = X(1:obj.k_delta,:); X_b2 = X(obj.k_delta+1:end,:);
                delta = theta(1:obj.k_delta); beta2 = theta(obj.k_delta+1:obj.k_delta+obj.k_lambda);
                pi = theta(obj.k_delta+obj.k_lambda+1:end);
                %
                X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
                for i = 1:obj.k_lambda
                    X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
                end
                %
                temp1 = obj.fcn_dg(X(obj.k_delta,:), pi(1));
                temp2 = X_b2;
                for i = 1:obj.k_lambda
                    temp2(i,:) = obj.fcn_dg(X(obj.k_delta+i,:), pi(i+1));
                end
                %
                if type_exclude_beta == 0
                    temp1 = delta(end) * temp1;
                    temp2 = repmat(beta2, 1, obj.n) .* temp2;
                elseif type_exclude_beta == 1
                    %do nothing
                end
                D_nu_hat = -[X_d; X_b2; temp1; temp2];
            end
        end
% %         function D_nu_hat = fcn_D_nu_hat_full(obj, theta, X, type_exclude_beta)
% %             num_var = size(X,1);
% %             X_d = X(1:obj.k_delta,:); 
% %             delta = theta(1:obj.k_delta); 
% %             if num_var == obj.k_delta  % only doing this for the full model b/c the null imposed model is only estimated with the full model estimation procedure
% %                 pi = theta(obj.k_delta+1:end);
% %                 X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
% %                 temp1 = obj.fcn_dg(X(obj.k_delta,:), pi(1));
% %                 if type_exclude_beta == 0
% %                     temp1 = delta(end) * temp1;
% %                 elseif type_exclude_beta == 1 %do nothing
% %                 end
% %                 D_nu_hat = -[X_d; temp1];
% %             elseif num_var > obj.k_delta
% %                 X_b2 = X(obj.k_delta+1:end,:);
% %                 beta2 = theta(obj.k_delta+1:obj.k_delta+obj.k_lambda);
% %                 pi = theta(obj.k_delta+obj.k_lambda+1:end);
% %                 %
% %                 X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
% %                 for i = 1:obj.k_lambda
% %                     X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
% %                 end
% %                 %
% %                 temp1 = obj.fcn_dg(X(obj.k_delta,:), pi(1));
% %                 temp2 = X_b2;
% %                 for i = 1:obj.k_lambda
% %                     temp2(i,:) = obj.fcn_dg(X(obj.k_delta+i,:), pi(i+1));
% %                 end
% %                 %
% %                 if type_exclude_beta == 0
% %                     temp1 = delta(end) * temp1;
% %                     temp2 = repmat(beta2, 1, obj.n) .* temp2;
% %                 elseif type_exclude_beta == 1
% %                     %do nothing
% %                 end
% %                 D_nu_hat = -[X_d; X_b2; temp1; temp2];
% %             end
% %         end
        function D_nu_hat = fcn_D_nu_hat_pars(obj, theta, X, type_exclude_beta, Taylor_flag)
            X_d = X(1:obj.k_delta,:); X_b2 = X(obj.k_delta+1,:);
            delta = theta(1:obj.k_delta); beta2 = theta(obj.k_delta+1);
            pi = theta(obj.k_delta+2:end);
            X_d(end,:) = obj.fcn_g(X_d(end,:), pi(1));            
            i = 1;
            X_b2(i,:) = obj.fcn_g(X_b2(i,:), pi(i+1));
            temp1 = obj.fcn_dg(X(obj.k_delta,:), pi(1));
            temp2 = obj.fcn_dg(X(obj.k_delta+1,:), pi(2));
            if type_exclude_beta == 0
                temp1 = delta(end) * temp1;
                temp2 = beta2 * temp2;
            elseif type_exclude_beta == 1
                %do nothing
            end
            D_nu_hat = -[X_d; X_b2; temp1; temp2];
        end
        %
        function m = fcn_m(obj, theta, X, Y, type)
            nu_hat = obj.fcn_nu_hat(theta, X, Y, type);
            m = nu_hat.^2 / 2;
        end
%         function Dm = fcn_Dm(obj, theta, X, Y, type, type_exclude_beta)
%             if nargin < 6
%                 type_exclude_beta = 0;
%             end
%             nu_hat = obj.fcn_nu_hat(theta, X, Y, type);
%             D_nu_hat = obj.fcn_D_nu_hat(theta, X, type, type_exclude_beta);
%             Dm = repmat(nu_hat,size(D_nu_hat,1),1) .* D_nu_hat;
%         end
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
                l_K_beta = l_K - (obj.k_delta + obj.k_lambda);
            elseif type == 1 % Pars model so only 1 lambda
                l_K_beta = l_K - (obj.k_delta + 1);
            end
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
            % chi = - (1/2) * tau' * inv(H) * tau;  % oops!
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
            if type == 0
                X = [obj.x_delta; obj.x_beta2];
%                 theta = obj.theta_hat_full;        
                theta = obj.theta_0;
                theta_0_in = obj.theta_0;
                i = [];
                l_K = obj.fcn_l_K_selection(theta, X, Y, type, i, type_exclude_beta);
                [mathfrac_z, HG]= obj.fcn_mathfrac_z(theta, X, Y, type, l_K, type_exclude_beta, theta_0_in, Z);
                V = HG * HG';
            elseif type == 1
                mathfrac_z = zeros((obj.k_delta+1+2)*obj.k_lambda,1);
                HG = zeros((obj.k_delta+1+2)*obj.k_lambda, obj.n);
                for i = 1:obj.k_lambda
                    X = [obj.x_delta; obj.x_beta2(i,:)];
                    ind_theta_i = [(1:obj.k_delta) obj.k_delta+i obj.k_delta+obj.k_lambda+1 obj.k_delta+obj.k_lambda+1+i];
%                     theta = obj.theta_hat_pars(:,i);        
                    theta = obj.theta_0(ind_theta_i);
                    theta_0_in = obj.theta_0(ind_theta_i);
                    l_K = obj.fcn_l_K_selection(theta, X, Y, type, i, type_exclude_beta);
                    ind = ( (obj.k_delta+1+2)*(i-1)+1:(obj.k_delta+1+2)*(i) );
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
            lambda_ind_full = (obj.ind_lambda_full);
            obj.lambda_distr_full = obj.mathfrac_z_distr_full(lambda_ind_full,:);
            lambda_ind_pars = ( obj.ind_lambda_pars : (obj.k_delta+1+1+1) : (obj.k_delta+1+2)*(obj.k_lambda) );
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
                ind_beta = (obj.k_delta:obj.k_delta+obj.k_lambda);
                ind_pi = ind_beta + (obj.k_lambda + 1);
                theta_temp = obj.theta_hat_full;
            elseif type == 1 % pars
                ind_beta = [obj.k_delta obj.k_delta+1];
                ind_pi = ind_beta + (1 + 1);
                theta_temp = obj.theta_hat_pars(:,i);
            end
            l_K_in = [];
            Z = ones(1, obj.n);
            G = obj.fcn_G(theta_temp, X, Y, type, l_K_in, type_exclude_beta, Z);
            H = obj.fcn_H(theta_temp, X, Y, type, l_K_in, type_exclude_beta);
            Sigma = inv(H) * (G * G') * inv(H)';
%             inv_Sigma = inv(Sigma);
            l_K_out = [];
            for j = 1:length(ind_beta) % for this experiment, H_0: pi = lambda = 0, so we test all the betas for inclusion in l_K
                % this differs in other experiments
                beta_j = theta_temp(ind_beta(j));
                ICS = sqrt( obj.n * beta_j * ( 1/Sigma(j,j) ) * beta_j );
                if (ICS <= obj.kappa_n) == 1
                    l_K_out = [l_K_out ind_pi(j)]; % include j's associate pi
                elseif (ICS > obj.kappa_n) == 1
                    % exclude j: do nothing here
                end
            end
            % l_K_out = [l_K_out ind_pi(2:end)]; % when H_0: beta = lambda = 0, we put all the lambdas in l_K
        end
        %
        % Estimation Function below
        function theta_out = fcn_estimation(obj, type, X_in, Y_in, bounds)
            if type == 0
                theta_out = obj.fcn_estimation_full(type, X_in, Y_in, bounds);
            elseif type == 1
                theta_out = obj.fcn_estimation_pars(type, X_in, Y_in, bounds);
            end
        end
        % Estimation via fmincon here
        function [Q, DQ] = fcn_loss(obj, theta, X, Y, type)
            m = obj.fcn_m(theta, X, Y, type);
            Dm = obj.fcn_Dm(theta, X, Y, type);
            Q = sum(m) / obj.n;
            DQ = sum(Dm,2) / obj.n;
        end
        function theta_out = fcn_estimation_full_2(obj, type, X_in, Y_in, bounds_pi)
            theta_init = [obj.delta_0; obj.beta2_0; obj.pi_0];
            LB = [repmat(-10,obj.k_delta+obj.k_lambda,1); bounds_pi(:,1)]; UB = [repmat(10,obj.k_delta+obj.k_lambda,1); bounds_pi(:,2)];
            [theta_out, fval_temp, exit_flag_temp] = fmincon(@(theta_in) obj.fcn_loss(theta_in, X_in, Y_in, type), theta_init, [], [], [], [], LB, UB, [], obj.options_fmincon);
        end
        function theta_out = fcn_estimation_pars_2(obj, type, X_in, Y_in, bounds_pi)
            theta_init = [obj.delta_0; obj.beta2_0; obj.pi_0];
            LB_all = [repmat(-10,obj.k_delta+obj.k_lambda,1); bounds_pi(:,1)]; UB_all = [repmat(10,obj.k_delta+obj.k_lambda,1); bounds_pi(:,2)];
            for i = 1:obj.k_lambda
                ind_pi = [1 1+i]; 
                ind_pars = [1:obj.k_delta obj.k_delta+i obj.k_delta+obj.k_lambda+ind_pi];
                theta_init = [obj.delta_0; obj.beta2_0(i); obj.pi_0(ind_pi)];
                LB = LB_all(ind_pars); UB = UB_all(ind_pars);
                X_pars = [X_in(1:obj.k_delta,:); X_in(obj.k_delta+i,:)];
                [theta_out, fval_temp, exit_flag_temp] = fmincon(@(theta_in) obj.fcn_loss(theta_in, X_pars, Y_in, type), theta_init, [], [], [], [], LB, UB, [], obj.options_fmincon);
            end
        end
        % Estimation via Levenberg-Marquardt Algorithm below
        function theta_out = fcn_estimation_full(obj, type, X_in, Y_in, bounds_pi)
            LB = bounds_pi(:,1); UB = bounds_pi(:,2);
            num_var = size(X_in,1);
            if num_var == obj.k_delta
                theta_init = [obj.delta_0; obj.pi_0(1)];
                R = zeros(1,obj.k_delta+1);
                    R(1, obj.k_delta+1) = 1;
            elseif num_var > obj.k_delta
                theta_init = [obj.delta_0; obj.beta2_0; obj.pi_0];
                R = zeros(length(obj.pi_0),length(obj.theta_0));
                for k = 1:length(obj.pi_0)
                    R(k, obj.k_delta+obj.k_lambda+k) = 1;
                end
            end
            % R = [0 0 1 0; 0 0 0 1]; 
            R = [R; -R];
            c = [LB; -UB]; % c = [-1; -1; -1; -1]; % LB, -UB
            theta_out = obj.fcn_estimation_proc(theta_init, X_in, Y_in, type, R, c);
        end
        function theta_out = fcn_estimation_pars(obj, type, X_in, Y_in, bounds_pi)
            theta_out = zeros(obj.k_delta+1+2, obj.k_lambda);
            LB = bounds_pi(:,1); UB = bounds_pi(:,2);
            for i = 1:obj.k_lambda
                % i=1; 
                ind_pi = [1 1+i]; 
                % ind_pi = [(1:obj.k_delta) obj.k_delta + i];
                theta_init = [obj.delta_0; obj.beta2_0(i); obj.pi_0(ind_pi)];
                ind_R = (obj.k_delta + 2:length(theta_init)); 
                % X_pars = [obj.x_delta(end,:); obj.x_beta2(i,:)];
                X_pars = [X_in(1:obj.k_delta,:); X_in(obj.k_delta+i,:)];
                R = zeros(length(ind_R),length(theta_init));
                for k = 1:length(ind_R)
                    R(k, ind_R(k)) = 1;
                end
                % R = [0 0 1 0; 0 0 0 1]; 
                R = [R; -R];
                c = [LB(ind_pi); -UB(ind_pi)]; % c = [-1; -1; -1; -1]; % LB, -UB
            %     theta_init = [.5; .0; 0; 0];
                theta_out(:,i) = obj.fcn_estimation_proc(theta_init, X_pars, Y_in, type, R, c);
            end
        end
        function theta_out = fcn_estimation_proc(obj, theta_init, X, Y, type, R, c)
%             type = 0; 
            tol = .00000001; check = 10000; it = 0;
            theta_in = theta_init; lambda_in = 1;
            f_in = obj.fcn_nu_hat(theta_in, X, Y, type);
            Df_in = obj.fcn_D_nu_hat(theta_in, X, type);
            while (it < 20) % (check > tol && it < 20)
                [theta_out, f_out, Df_out, lambda_out] = obj.fcn_NLLS(theta_in, lambda_in, Df_in, f_in, X, Y, type, R, c);
                check = abs(theta_in - theta_out)' * abs(theta_in - theta_out);
                it = it + 1;
                theta_in = theta_out; f_in = f_out; Df_in = Df_out; lambda_in = lambda_out;
            end
        end
        function [theta_out, f_out, Df_out, lambda_out] = fcn_NLLS(obj, theta_in, lambda_in, Df_in, f_in, X, Y, type, R, c)
            % Levenberg-Marquardt Algorithm
            A = ( Df_in * Df_in' + (lambda_in .* eye(size(Df_in,1)) ) ) \ eye(size(Df_in,1));
            B = Df_in * f_in';
            theta_out = theta_in - A * B;
            %
            if nargin > 7
                mu = zeros(size(R,1), 1);
                % ind = find(R * theta_out - c > 0);
                ind2 = find(R * theta_out - c <= 0);
                % R_tilde = R; R_tilde(ind, :) = [];
                % c_tilde = c; c_tilde(ind, :) = [];
                R_tilde = R(ind2, :);
                c_tilde = c(ind2);
                D = ( R_tilde * A * R_tilde') \ eye(size(R_tilde,1));
                E = ( c_tilde - R_tilde * theta_in) + R_tilde * A * B;
                temp = D * E;
                mu(ind2) = temp;
                theta_out = theta_in - A * ( B - R' * mu );
            end
            %
            f_out = obj.fcn_nu_hat(theta_out, X, Y, type);
            Df_out = obj.fcn_D_nu_hat(theta_out, X, type);
            norm_f_out = f_out * f_out'; norm_f_in = f_in * f_in';
            if norm_f_out < norm_f_in
                lambda_out = .9 * lambda_in;
            elseif norm_f_out >= norm_f_in
                lambda_out = 2 * lambda_in;
                theta_out = theta_in; f_out = f_in; Df_out = Df_in;
            end
%             theta_in = theta_out; f_in = f_out; Df_in = Df_out;
        end
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

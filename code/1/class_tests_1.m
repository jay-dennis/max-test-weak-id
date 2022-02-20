classdef class_tests_1 < class_dgp_1
    properties (SetAccess=public, GetAccess=public)
        %% Max test and max t-test variables
        k_lambda_n
        theta_star
        y_hat
        % V_all, V_theta
        max_test_stat, max_t_test_stat, % max_t_test_stat_orig
        max_temp, max_t_temp
        % T_max_test, T_max_t_test
        p_max_test, p_max_t_test
        p_max_test_bs, p_max_t_test_bs
        dr_max
        dr_max_t
        M_max_sim
        %% Asymptotic Wald Test
        wald_Cheng_test_stat
%         wald_orig_test_stat
        T_wald
        p_wald
        %% Parametric Bootstrapped Wald Tests
        M
        % T_star_wald_bs
        p_wald_bs
        %% Taylor Expansion Variables
        theta_hat_full_Taylor, theta_hat_pars_Taylor
        lambda_hat_full_Taylor, lambda_hat_pars_Taylor
        Wald_Taylor, Max_Taylor, Max_t_Taylor
        theta_hat_full_Taylor_bs, theta_hat_pars_Taylor_bs
        lambda_hat_full_Taylor_bs, lambda_hat_pars_Taylor_bs
        Wald_Taylor_bs, Max_Taylor_bs, Max_t_Taylor_bs
        %% other
        alpha_levels
        dr
        p_vals
    end
    properties (SetAccess=public, GetAccess=public, Hidden = true)
    end
    methods
        function out = about(obj)
            fprintf('\n \n Generates \n');
            fprintf('Version 0 \n');
            fprintf('------------ \n');
        end
        function obj = class_tests_1(dgp_type, hypothesis_type, n, seed_input, k_lambda, k_delta, beta_in, k_lambda_n)
            obj@class_dgp_1(dgp_type, hypothesis_type, n, seed_input, k_lambda, k_delta, beta_in);
            obj.k_lambda_n = k_lambda_n;
        end
        function obj = set_params(obj, k_lambda_n)
            obj.alpha_levels = [0.01 0.05 0.10];
            obj.k_lambda_n = k_lambda_n;
            obj.M_max_sim=25000; % 25000;
            obj.M = 1000;
        end
        function obj = clean_up(obj)
            obj.Y0=[];
            obj.e=[];
            obj.x_delta=[];
            obj.x_lambda=[];
            obj.x_pi=[];
            obj.y_hat=[];
            obj.options_fmincon = [];
            obj.init_n = [];
            % obj.theta_0 = [];
            % obj.delta_0 = [];
        end
        % General Functions
        function dr_vec = set_decision_rules(obj)
            p_temp = obj.p_vals;
            len   = length(p_temp);
            len_a = length(obj.alpha_levels);
            dr_vec = zeros(len,len_a);
            for a = 1:len_a
                alpha_level = obj.alpha_levels(a);
                for l = 1:len
                    p = p_temp(l);
                    if p >= alpha_level
                        dr_temp = 0; % Fail to reject H_0
                    elseif p < alpha_level
                        dr_temp = 1; % Reject H_0
                    end
                    dr_vec(l,a) = dr_temp;
                end
            end
        end
        function p = pval_fcn(obj, T, M, statistic)
            p = (1/M) * sum((T>statistic));
        end
        %
        %%%%%%%%%%%%%
        %%%%%%%%%%%%%
        % Initial Tests
        function obj = fcn_tests_initial(obj)
            Y = obj.Y0;
            X = [obj.x_delta; obj.x_lambda];

            type_exclude_beta = 0;
            type = 0; theta_in = obj.theta_hat_full;
            [V_hat_full, V_hat_full_lambda] = obj.fcn_V_hat(theta_in, X, Y, type, type_exclude_beta);
            type = 1; theta_in = obj.theta_hat_pars;
            [V_hat_pars, V_hat_pars_lambda] = obj.fcn_V_hat(theta_in, X, Y, type, type_exclude_beta);
            
            % subtracting lambda_0 for initial observations only
%             obj.wald_Cheng_test_stat = obj.n * (obj.lambda_hat_full - obj.lambda_0)' * inv(V_hat_full_lambda) * (obj.lambda_hat_full - obj.lambda_0);
            obj.wald_Cheng_test_stat = obj.n * (obj.lambda_hat_full)' * inv(V_hat_full_lambda) * (obj.lambda_hat_full);

            W = ones(obj.k_lambda_n,1);
            Wt = sqrt(diag(V_hat_pars_lambda)).^(-1);
%             obj.max_test_stat = sqrt(obj.n) * max(abs(W .* (obj.lambda_hat_pars - obj.lambda_0)));
%             obj.max_t_test_stat = sqrt(obj.n) * max(abs(Wt .* (obj.lambda_hat_pars - obj.lambda_0)));
            obj.max_test_stat = sqrt(obj.n) * max(abs(W .* (obj.lambda_hat_pars)));
            obj.max_t_test_stat = sqrt(obj.n) * max(abs(Wt .* (obj.lambda_hat_pars)));
            
            % Simulate max test with correct covar matrix here for standard
            % critical values
                centering = zeros(size(obj.theta_hat_pars,1),size(obj.theta_hat_pars,2)); 
                centering(1,:) = obj.theta_hat_pars(1,:); % impose the null
                theta_sim = randn((obj.k_delta+1+2)*obj.k_lambda,1); 
                theta_sim = sqrtm(V_hat_pars) * theta_sim + centering(:);
                ind = ((obj.k_delta+1):(obj.k_delta+1+2):(obj.k_delta+1+2)*obj.k_lambda);
                lambda_sim = theta_sim(ind);
                obj.max_temp = max(abs(lambda_sim))';
                obj.max_t_temp = max(abs(Wt .* lambda_sim))';
            
                
            obj = obj.fcn_linearized_tests;
                
               
%                 lambda_sim = randn(obj.k_lambda,1); 
%                 obj.max_t_temp = max(abs(lambda_sim))';
                 
%             type_exclude_beta = 0;
%             [~, V_hat_full_lambda_orig] = obj.fcn_V_hat(X, Y, 0, type_exclude_beta);
%             [~, V_hat_pars_lambda_orig] = obj.fcn_V_hat(X, Y, 1, type_exclude_beta);
%             obj.wald_orig_test_stat = obj.n * obj.lambda_hat_full' * inv(V_hat_full_lambda_orig) * obj.lambda_hat_full;
%             Wt_orig = sqrt(diag(V_hat_pars_lambda_orig)).^(-1);
%             obj.max_t_test_stat_orig = sqrt(obj.n) * max(abs(Wt_orig .* obj.lambda_hat_pars));   
        end
        function [V_hat, V_hat_lambda] = fcn_V_hat(obj, theta_in, X_in, Y, type, type_exclude_beta)
            l_K = []; % use ICS pre-test to determine this for the critical values
            if type == 0
                G_hat_full = obj.fcn_G(theta_in, X_in, Y, type, l_K, type_exclude_beta);
                H_hat_full = obj.fcn_H(theta_in, X_in, Y, type, l_K, type_exclude_beta);
                V_hat = inv(H_hat_full) * (G_hat_full * G_hat_full') * inv(H_hat_full)';
                lambda_ind = [(obj.k_delta + 1) : (obj.k_delta + obj.k_lambda)];
                V_hat_lambda = V_hat(lambda_ind, lambda_ind);
            elseif type == 1
                A = zeros(obj.k_lambda * (obj.k_delta+1+2), obj.n);
                lambda_ind = zeros(obj.k_lambda, 1);
                for i = 1:obj.k_lambda
                    X = [X_in(1:obj.k_delta,:); X_in(obj.k_delta+i,:)];
                    l_K = [];
                    G_i = obj.fcn_G(theta_in(:,i), X, Y, type, l_K, type_exclude_beta); 
                    H_i = obj.fcn_H(theta_in(:,i), X, Y, type, l_K, type_exclude_beta); 
                    A_i = inv(H_i) * G_i;
                    ind = ( ((i-1) * (obj.k_delta+1+2)+1) : (i*(obj.k_delta+1+2)) );
                    A(ind,:) = A_i;
                    lambda_ind(i) = (i-1) * (obj.k_delta+1+2)+2;
                end
                V_hat = A * A';
                V_hat_lambda = V_hat(lambda_ind, lambda_ind);
            end
        end
        %
        function obj = fcn_linearized_tests(obj)
            [obj.theta_hat_full_Taylor, obj.theta_hat_pars_Taylor, V_full, V_pars] = obj.fcn_linearized_estimation;
            
            len_d = obj.k_delta+1;
            len_l = 2*obj.k_lambda;
            obj.lambda_hat_full_Taylor = obj.theta_hat_full_Taylor(len_d+1:end);
            obj.lambda_hat_pars_Taylor = obj.theta_hat_pars_Taylor(end,:)';
            V_full_lambda = V_full(len_d+1:len_d+len_l, len_d+1:len_d+len_l);
            ind = ( len_d+1 : len_d+1 : (len_d+1)*len_l );
            V_pars_lambda = V_pars(ind, ind);
            
            obj.Wald_Taylor = obj.n * obj.lambda_hat_full_Taylor' * inv(V_full_lambda) * obj.lambda_hat_full_Taylor;
            obj.Max_Taylor = sqrt(obj.n) * max( abs( obj.lambda_hat_pars_Taylor ) )';
            Wt = sqrt( diag(V_pars_lambda) ).^(-1);
            obj.Max_t_Taylor = sqrt(obj.n) * max( abs( Wt .* obj.lambda_hat_pars_Taylor ) )';
        end
        function [theta_hat_full_Taylor, theta_hat_pars_Taylor, V_full, V_pars] = fcn_linearized_estimation(obj)
            % Estimates a Taylor Approximation to the nonlinear model
            % Here we Linearize every nonlinear component of the model,
            % including those associated with parameters that will not be
            % tested.  This is for convenience only.
            
            % First construct the regressors
            Y = obj.Y0;
            
            X_d = obj.x_delta; 
            X_d2 = obj.x_delta;
                X_d(end,:) = obj.fcn_g(obj.x_delta(end,:), 0);
                X_d2(end,:) = obj.fcn_dg(obj.x_delta(end,:), 0);
            
            % Construct the new regressors from a Taylor expansion 
            X_l = obj.x_lambda;
            X_l2 = obj.x_lambda;
            for i = 1:obj.k_lambda
                X_l(i,:) = obj.fcn_g(obj.x_lambda(i,:), 0);
                X_l2(i,:) = obj.fcn_dg(obj.x_lambda(i,:), 0);
            end
            
            X_d = [X_d; X_d2];          
            X_l = [X_l; X_l2];
            X_all = [X_d; X_l];
            
            % Full
            theta_hat_full_Taylor = (Y / X_all)';
            % Covar matrix
            e_hat_full = Y - (theta_hat_full_Taylor' * X_all);
            H = X_all * X_all' / obj.n; 
            G = X_all .* repmat(e_hat_full, size(X_all,1), 1) / sqrt(obj.n);
            V_full = inv(H) * (G * G') * inv(H);
            
            % Pars
            len_d = size(X_d,1); len_l = size(X_l,1);
            theta_hat_pars_Taylor = zeros( len_d+1, len_l );
            temp = zeros((len_d+1)*len_l, obj.n);
            for i = 1:len_l
                X_i = [X_d; X_l(i,:)];
                temp_theta = (Y / X_i)';
                theta_hat_pars_Taylor(:,i) = temp_theta;
                e_hat_i = Y - (temp_theta' * X_i);
                Hi = X_i * X_i' / obj.n; 
                Gi = X_i .* repmat(e_hat_i, len_d+1, 1) / sqrt(obj.n);
                V_i = inv(Hi) * Gi ;
                ind = (i-1)*(len_d+1)+1:(i)*(len_d+1);
                temp(ind,:) = V_i;
            end
            V_pars = temp * temp';
        end
        
        function obj = fcn_linearized_bs(obj, M, seed)
            obj.cluster_bs2_distr_seed = seed;
            rng(seed);
            Z_mat = randn(M, obj.n);

            Y = obj.Y0;
            X_d = obj.x_delta; 
            X_d2 = obj.x_delta;
                X_d(end,:) = obj.fcn_g(obj.x_delta(end,:), 0);
                X_d2(end,:) = obj.fcn_dg(obj.x_delta(end,:), 0);
            
            % Construct the new regressors from a Taylor expansion 
            X_l = obj.x_lambda;
            X_l2 = obj.x_lambda;
            for i = 1:obj.k_lambda
                X_l(i,:) = obj.fcn_g(obj.x_lambda(i,:), 0);
                X_l2(i,:) = obj.fcn_dg(obj.x_lambda(i,:), 0);
            end
            
            X_d = [X_d; X_d2];          
            X_l = [X_l; X_l2];
            X_all = [X_d; X_l];
            len_d = size(X_d,1); len_l = size(X_l,1);

            Wald_Taylor_bs_temp = zeros(M,1);
            Max_Taylor_bs_temp = zeros(M,1);
            Max_t_Taylor_bs_temp = zeros(M,1);
            obj.theta_hat_full_Taylor_bs  = zeros(len_d+len_l,M);
            obj.theta_hat_pars_Taylor_bs  = zeros((len_d+1)*len_l,M);
            obj.lambda_hat_full_Taylor_bs = zeros(len_l,M);
            obj.lambda_hat_pars_Taylor_bs = zeros(len_l,M);
           
            %estimate the null imposed model
            theta_hat_Taylor_null = (Y / X_d)';
            e_hat_full = Y - (theta_hat_Taylor_null' * X_d);
            
            H = X_all * X_all' / obj.n; 
            for m = 1:M
                Z_m = Z_mat(m, :);
                Y_m = (theta_hat_Taylor_null' * X_d) + (e_hat_full .* Z_m);
    
                % Full
                theta_hat_full_Taylor_temp = (Y_m / X_all)';
                % Covar matrix
                e_hat_full = Y_m - (theta_hat_full_Taylor_temp' * X_all);
                G = X_all .* repmat(e_hat_full, size(X_all,1), 1) / sqrt(obj.n);
                V_full = inv(H) * (G * G') * inv(H);

                % Pars
                theta_hat_pars_Taylor_temp = zeros( len_d+1, len_l );
                temp = zeros((len_d+1)*len_l, obj.n);

                for i = 1:len_l
                    X_i = [X_d; X_l(i,:)];
                    Hi = X_i * X_i' / obj.n; 
                    
                    temp_theta = (Y_m / X_i)';
                    theta_hat_pars_Taylor_temp(:,i) = temp_theta;
                    e_hat_i = Y_m - (temp_theta' * X_i);
                    Gi = X_i .* repmat(e_hat_i, len_d+1, 1) / sqrt(obj.n);
                    V_i = inv(Hi) * Gi ;
                    ind = (i-1)*(len_d+1)+1:(i)*(len_d+1);
                    temp(ind,:) = V_i;
                end
                V_pars = temp * temp';
                
                %construct test stat
                lambda_hat_full_Taylor_temp = theta_hat_full_Taylor_temp(len_d+1:end);
                V_full_lambda = V_full(len_d+1:len_d+len_l, len_d+1:len_d+len_l);
            
                lambda_hat_pars_Taylor_temp = theta_hat_pars_Taylor_temp(end,:)';
                ind = ( len_d+1 : len_d+1 : (len_d+1)*len_l );
                V_pars_lambda = V_pars(ind, ind);
   
                Wald_Taylor_temp = obj.n * lambda_hat_full_Taylor_temp' * inv(V_full_lambda) * lambda_hat_full_Taylor_temp;
                Max_Taylor_temp = sqrt(obj.n) * max( abs( lambda_hat_pars_Taylor_temp ) )';
                Wt = sqrt( diag(V_pars_lambda) ).^(-1);
                Max_t_Taylor_temp = sqrt(obj.n) * max( abs( Wt .* lambda_hat_pars_Taylor_temp ) )';
                
                % store the test stats
                Wald_Taylor_bs_temp(m) = Wald_Taylor_temp;
                Max_Taylor_bs_temp(m) = Max_Taylor_temp;
                Max_t_Taylor_bs_temp(m) = Max_t_Taylor_temp;

                obj.theta_hat_full_Taylor_bs(:,m)  = theta_hat_full_Taylor_temp;
                obj.theta_hat_pars_Taylor_bs(:,m)  = theta_hat_pars_Taylor_temp(:);
                obj.lambda_hat_full_Taylor_bs(:,m) = lambda_hat_full_Taylor_temp;
                obj.lambda_hat_pars_Taylor_bs(:,m) = lambda_hat_pars_Taylor_temp;
            end
            
        
            obj.Wald_Taylor_bs = Wald_Taylor_bs_temp;
            obj.Max_Taylor_bs = Max_Taylor_bs_temp;
            obj.Max_t_Taylor_bs = Max_t_Taylor_bs_temp;

        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %%%%%%%%%%%%%
        %%%%%%%%%%%%%
        %%%%%%%%%%%%%
        % All tests
        function obj = run_all_tests_fcn(obj, k_lambda_n)
            % k_lambda_n = 10; %remove this later
            obj = obj.set_params(k_lambda_n);
            obj = obj.max_test_stat_fcn;
            %
            [obj.T_wald, obj.p_wald] = obj.wald_fcn;
            [~, obj.p_wald_bs] = obj.bs_wald_fcn(obj.T_wald);
            %
            obj.p_vals = [obj.p_max_test; obj.p_max_t_test; obj.p_max_test_bs; obj.p_max_t_test_bs; obj.p_wald; obj.p_wald_bs];
            obj.dr = obj.set_decision_rules;
        end
        %
        %%%%%%%%%%%%%
        % Max Test Methods
        function obj = max_test_stat_fcn(obj)
            obj.theta_star = zeros(obj.k_lambda_n,(obj.k_delta + 1));
            obj.y_hat = zeros(obj.n,obj.k_lambda_n);
            for i = 1:obj.k_lambda_n
                temp_x = [obj.x_delta; obj.x_theta(i, :)];
                temp_theta = temp_x' \ obj.Y0'; % Just OLS here;
                obj.theta_star(i,:) = temp_theta;
                obj.y_hat(:,i) = temp_theta' * temp_x;
            end
            [V_all, V_theta] = obj.V_hat_fcn;
            
            W = ones(obj.k_lambda_n,1);
            obj.max_test_stat = sqrt(obj.n) * max(abs(W .* obj.theta_star(:,end)));
            Wt = sqrt(diag(V_theta)).^(-1);
            obj.max_t_test_stat = sqrt(obj.n) * max(abs(Wt .* obj.theta_star(:,end)));
            
            W_in = [W Wt]; stat_in = [obj.max_test_stat obj.max_t_test_stat];
            [~, p_out] = obj.max_test_pval_fcn(V_all, obj.M_max_sim, W_in, stat_in);
            obj.p_max_test = p_out(1); obj.p_max_t_test = p_out(2);
            
            [~, p_out] = obj.max_test_bs_pval_fcn(obj.M_max_sim, W_in, stat_in);
            obj.p_max_test_bs = p_out(1); obj.p_max_t_test_bs = p_out(2);
        end
        function [V, V_theta] = V_hat_fcn(obj)
            V = zeros(obj.k_lambda_n*(obj.k_delta + 1), obj.k_lambda_n*(obj.k_delta + 1));
            V_theta = zeros(obj.k_lambda_n,obj.k_lambda_n);
            for i = 1:obj.k_lambda_n
                x_i = [obj.x_delta; obj.x_theta(i, :)];
                delta_hat_i = obj.theta_star(i,1:end-1)';
                J_i = (x_i * x_i') / obj.n;
                temp_e_hat_i = obj.Y0 - delta_hat_i' * obj.x_delta;
                for j = 1:obj.k_lambda_n
                    x_j = [obj.x_delta; obj.x_theta(j, :)];
                    delta_hat_j = obj.theta_star(j,1:end-1)';
                    J_j = (x_j * x_j') / obj.n;
                    temp_e_hat_j = obj.Y0 - delta_hat_j' * obj.x_delta;
                    temp_i = repmat(temp_e_hat_i, obj.k_delta+1, 1) .* x_i;
                    temp_j = repmat(temp_e_hat_j, obj.k_delta+1, 1) .* x_j;
                    S = temp_i * temp_j' / obj.n; 
                    v = ( J_i \ eye(size(J_i,1)) ) * S * ( J_j \ eye(size(J_j,1)) );
                    V(((i-1)*(obj.k_delta + 1)+1):(i*(obj.k_delta + 1)),((j-1)*(obj.k_delta + 1)+1):(j*(obj.k_delta + 1))) = v;
                    V_theta(i,j) = v(end,end);
                end
            end
        end
        function [T_out, p_out] = max_test_pval_fcn(obj, V, M, W_in, stat_in)
            %temp_v = chol(V);
            temp_v = sqrtm(V);
            Z = randn(obj.k_lambda_n*(obj.k_delta + 1),M);
            Z_sim = temp_v * Z;
            Z_theta = zeros(obj.k_lambda_n,M);
            for i = 1:obj.k_lambda_n
                Z_theta(i,:) = Z_sim(i*(obj.k_delta + 1),:);
            end
            T_out = zeros(length(stat_in), M); p_out = zeros(length(stat_in),1);
            for j = 1:length(stat_in)
                W = W_in(:, j); statistic = stat_in(j);
                W = repmat(W, 1, M);
                T = max(abs(W .* Z_theta),[],1);
                p = obj.pval_fcn(T, M, statistic);
                T_out(j,:) = T; p_out(j) = p;
            end
        end
        function [T_out, p_out] = max_test_bs_pval_fcn(obj, M, W_in, stat_in)
            M = 1000;
            % Bootstrapped Max Test
            % estimate imposing the null.
            temp_X = obj.x_delta';
            % d_hat_null = inv(temp_X' * temp_X) * temp_X' * obj.Y0';
            d_hat_null = temp_X \ obj.Y0'; % faster
            % form null imposed residuals.
            y_hat_null = d_hat_null' * obj.x_delta;
            e_hat_null = obj.Y0 - y_hat_null;
            % mutliply resids by normal rv
            Z = randn(M, obj.n);
            e_hat_wbs = repmat(e_hat_null, M, 1) .* Z;
            % construct y*
            y_wbs = repmat(y_hat_null, M, 1) + e_hat_wbs;
            % estimate pars models
            theta_hat_wbs = zeros(M, obj.k_lambda_n);
            for m = 1:M
                temp_Y = y_wbs(m,:)';
                for i = 1:obj.k_lambda_n
                    temp_X = [obj.x_delta' obj.x_theta(i,:)'];
                    % b_hat_wbs = inv(temp_X' * temp_X) * temp_X' * temp_Y;
                    b_hat_wbs = temp_X \ temp_Y;
                    theta_hat_wbs(m,i) = b_hat_wbs(end);
                end
            end
            % form max stats, and sort
            % feed to pval func
            T_out = zeros(M, length(stat_in)); p_out = zeros(length(stat_in),1);
            for j = 1:length(stat_in)
                W = W_in(:, j); statistic = stat_in(j);
                W = repmat(W', M, 1);
                T = sqrt(obj.n) * max(abs(W .* theta_hat_wbs),[],2);
                T = sort(T);
                p = obj.pval_fcn(T, M, statistic);
                T_out(:,j) = T; p_out(j) = p;
            end
        end
        %
        % Wald
        function [wald, p_val] = wald_fcn(obj)
            y0 = obj.Y0';
            X = [obj.x_delta' obj.x_theta(1:obj.k_lambda_n,:)'];
            temp_beta = X \ y0; % Just OLS here;
            e_hat = y0 -  X * temp_beta;
            % Use White's Robust variance estimator since e may contain
            % omitted variables.
            temp = X .* repmat(e_hat, 1, obj.k_delta + obj.k_lambda_n);
            % var_B_White = inv(X' * X) * (Z' * Z) * inv(X' * X);
            outer = ((X' * X) \ eye(obj.k_delta + obj.k_lambda_n));
            inner = temp' * temp;
            v = outer * inner * outer;
            v_theta = v((obj.k_delta+1):end,(obj.k_delta+1):end);
            v_inv = (v_theta \ eye(obj.k_lambda_n));
            %
            theta_vec = temp_beta((obj.k_delta+1):end);
            wald = (theta_vec)' * v_inv * (theta_vec);
            p_val = 1 - chi2cdf(wald, obj.k_lambda_n);
        end
        % Bootstrapped Wald
        function [T_star, p_val] = bs_wald_fcn(obj, wald)
            % impose the null first
            X_null = obj.x_delta'; X_all = [obj.x_delta' obj.x_theta(1:obj.k_lambda_n,:)'];
            y0 = obj.Y0';
            beta_hat_null = X_null \ y0; % Just OLS here;
            y_hat_null = X_null * beta_hat_null;
            e_null = y0 - y_hat_null;
            % generate null imposed y's
            Z = randn(obj.n, obj.M);
            e_bs = repmat(e_null, 1, obj.M) .* Z;
            y_null_bs = y_hat_null + e_bs;
            %
            T_star = zeros(obj.M,1);
            for m = 1:obj.M
                y_temp = y_null_bs(:,m);
                b_hat = X_all \ y_temp;
                theta_vec = b_hat(obj.k_delta+1:end);
                e_hat = y_temp - X_all * b_hat;
                %
                temp = X_all .* repmat(e_hat, 1, obj.k_delta + obj.k_lambda_n);
                outer = ((X_all' * X_all) \ eye(obj.k_delta + obj.k_lambda_n));
                inner = temp' * temp;
                v = outer * inner * outer;
                %
                v_theta = v((obj.k_delta+1):end,(obj.k_delta+1):end);
                v_inv = (v_theta \ eye(obj.k_lambda_n));
                T_star(m) = (theta_vec)' * v_inv * (theta_vec);
            end
            T_star = sort(T_star);
            p_val = obj.pval_fcn(T_star, obj.M, wald);
        end
        %
    end
end

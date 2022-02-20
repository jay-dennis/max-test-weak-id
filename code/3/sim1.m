function sim1(dgp_type, n, J, seed_type)

    % sim_2_dgp1_hyp1_n250
    if nargin == 0
        clear all; clc;
        J = 10;            %number of Simulations
        dgp_type = 1; %1:3
        n = 100;
        seed_type = 1;
    end

    test_number = 1;

    if seed_type == 0  % all at once, random
        seed_input = 'shuffle';
        rng('shuffle');
        temp = rng;
        sim_number = temp.Seed;
        clear temp;
    elseif seed_type == 1  % all at once, reproducible
        seed_input0 = 0;
        sim_number = seed_input0;
    elseif seed_type == 2  % longleaf array, reproducible
        seed_input0 = str2num(getenv('SLURM_ARRAY_TASK_ID'));
        seed_input0 = seed_input0 - 1;
        sim_number = seed_input0;
    end

    if n == 1000
        k_theta_vec = [40, floor(n/2)];
        k_delta_vec = [1, 20];
        k_theta_n_vec = [20, floor(5*n^(1/2 - .0001))];
    else
        k_theta_vec = [10, 40, floor(n/2)];
        k_delta_vec = [1, 20];
        k_theta_n_vec = [1, 10, 20, floor(5*n^(1/2 - .0001))];
    end
    
    beta1_vec = unique([0 1/sqrt(n) 1]);
    beta2_vec = unique([0 1/sqrt(n) 2/sqrt(n) 5/sqrt(n) 10/sqrt(n) 1]);
%     hypothesis_type_vec = (1:length(beta2_vec));
%     hypothesis_type_vec(1) = 1; hypothesis_type_vec(end) = 3; 
        
    warning('off','all');

    for k_delta = k_delta_vec
        for k_theta = k_theta_vec
            for beta1_ind = 1:length(beta1_vec)
                for beta2_ind = 1:length(beta2_vec)
                    beta1 = beta1_vec(beta1_ind); b1_type = beta1_ind;
                    beta2 = beta2_vec(beta2_ind); hypothesis_type = beta2_ind;
                    beta_in = [beta1 beta2];

                    data0 = [];

                    tic;
                    for j = 1:J
                        if (seed_type == 1 || seed_type == 2)
                            seed_input = seed_input0 * 10^ceil(log(J)/log(10)) + j;
                        end
                        data0 = [data0 class_tests_1(dgp_type, hypothesis_type, n, seed_input, k_theta, k_delta, beta_in, k_lambda_n)];
                    end
                    time.dgp = toc / 60;

                    for k_theta_n = k_theta_n_vec
                        fprintf('DGP: %d, Hyp: %d, k_d: %d, k_t: %d, k_t_n: %d \n', dgp_type, hypothesis_type, k_delta, k_theta, k_theta_n)
                        if k_theta_n <= k_theta
                            data = data0;
                            tic;
                            for j = 1:J
                                data(j) = data(j).run_all_tests_fcn(k_theta_n);
                            end
                            time.est = toc / 60;

                            %clean up
                            for j = 1:J
                                data(j) = data(j).clean_up();
                            end

                            %output_main_dir = sprintf('G:/Simulation_data/max_test_many_zeros/10July/data_n%d', n);
                            output_main_dir = sprintf('./data_n%d', n);

                            outputdir = sprintf('%s/output_%d_dgp%d_hyp%d_b1%d_n%d_kd%d_kt%d_ktn%d', output_main_dir, test_number, dgp_type, hypothesis_type, b1_type, n, k_delta, k_theta, k_theta_n);
                            if exist(outputdir,'dir') == 0
                                mkdir(outputdir);
                            end
                            outputname=sprintf('%s/data_%d_dgp%d_hyp%d_b1%d_n%d_ktn%d_J%d_%d.mat', outputdir, test_number, dgp_type, hypothesis_type, b1_type, n, k_theta_n, J, sim_number);

                            count = 0;
                            while exist(outputname,'file') == 2
                                fprintf('Error: Name Taken; Iterating Name Forward \n')
                                count = count + 1;
                                outputname=sprintf('%s/data_%d_dgp%d_hyp%d_b1%d_n%d_ktn%d_J%d_%d_%d.mat', outputdir, test_number, dgp_type, hypothesis_type, b1_type, n, k_theta_n, J, sim_number, count);
                            end
                            % display(outputname);
                            save(outputname, 'data', 'time');

                            clear data;

                        end % if k_theta_n <= k_theta
                    end % k_theta_n loop
                    clear data0 time;
                end
            end
        end % k_theta loop
    end % k_delta loop

end % function
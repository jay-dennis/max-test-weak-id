clear all;
clc;

cluster_job_type = 'array';

sim_file_name = 'f0_cluster';
cons_file_name = 'f2_consolodate_wrapper';
results_file_name = 'f3_results';
test_number = 1;
dgp_type_vec = 2; % 1:3

dgp_type = 2;
n = 500;
k_delta = 1;
k_lambda = 50;
beta1 = 1;
beta2 = 0;
hypothesis_type = 3;
seed_type = 2;

T_vec = n; % 200; % [100 200 500 1000]

starting_sim_number = 1; % be careful here.
total_num_sims = 10000;
num_sims_per_array_job = 25;
array_batch_size = 20; % number of jobs the cluster runs at one time. (the rest go in the queue.)


num_array_jobs = ceil(total_num_sims / num_sims_per_array_job);
array_start = floor(starting_sim_number / num_sims_per_array_job) + 1;
% array_len{100} = num_array_jobs; array_len{200} = num_array_jobs; array_len{250} = num_array_jobs; array_len{500} = num_array_jobs; array_len{1000} = num_array_jobs;
array_len = num_array_jobs; 

infile = fullfile(sprintf('%s.m', sim_file_name));

cluster_dir = sprintf('./cluster_sub');
if exist(cluster_dir,'dir') == 0
    mkdir(cluster_dir);
end

for n = T_vec
    if strcmp(cluster_job_type, 'array') == 1
        J = floor(total_num_sims / array_len);
        array_end = array_start + array_len - 1;
    end
    outputname_ll = sprintf('ll_%d_T%d', test_number, n);
    fname = sprintf('%s/%s.txt', cluster_dir, outputname_ll);
    FID = fopen(fname, 'w');
    %
    outputname_consolodate = sprintf('consolodate_%d_T%d', test_number, n);
    fname_cons = sprintf('%s/%s.txt', cluster_dir, outputname_consolodate);
    FID_cons = fopen(fname_cons, 'w');
    %
    outputname_results = sprintf('results_%d_T%d', test_number, n);
    fname_results = sprintf('%s/%s.txt', cluster_dir, outputname_results);
    FID_results = fopen(fname_results, 'w');
    %
    %allocate 10 min/sim (for longleaf) * J * T/100
    wall_time = 20*J*n/100;
    for dgp_type = dgp_type_vec
        job_name = sprintf('MTw%dd%dT%dJ%d', test_number, dgp_type, n, J);
        if strcmp(cluster_job_type, 'array') == 1
            m_file_name = sprintf('%s(%d,%d,%d,%d,%d,%d,%d,%d,%d)', sim_file_name, dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, seed_type);
            fprintf(FID, 'sbatch --job-name=%s -o %s.%%j.txt -t %d --array=[%d-%d%%%d] --wrap=\"matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r ''%s; quit;''\"; \n \n', job_name, job_name, wall_time, array_start, array_end, array_batch_size, m_file_name);
            %
            cons_file_name = sprintf('%s(%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)', cons_file_name, dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type, J, array_start, array_end);
            fprintf(FID_cons, '%s \n \n', cons_file_name);
            %
            results_file_name = sprintf('%s(%d,%d,%d,%d,%d,%d,%d)', results_file_name, dgp_type, n, k_delta, k_lambda, beta1, beta2, hypothesis_type);
            fprintf(FID_results, '%s \n \n', results_file_name);
        end
        fprintf(FID, '\n');
        fprintf(FID_cons, '\n');
        fprintf(FID_results, '\n');
    end
    fclose(FID);
    fclose(FID_cons);
    fclose(FID_results);
end

outfile = fullfile(cluster_dir, sprintf('%s.m', sim_file_name));
copyfile(infile, outfile);

tempfile = sprintf('class_dgp_%d.m', test_number);
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);

tempfile = sprintf('class_tests_%d.m', test_number);
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);

tempfile = sprintf('f0_1_sim_distr.m');
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);

tempfile = sprintf('f1_1_sims.m');
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);



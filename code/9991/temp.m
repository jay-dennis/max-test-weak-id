    lambda_hat_full = zeros(J,1);
    lambda_hat_pars = zeros(J,1);

    for j = 1:J
        lambda_hat_full(j) = data0(j).lambda_hat_full;
        lambda_hat_pars(j) = data0(j).lambda_hat_pars;
    end

    
    fig = fig+1;
    figure(fig);
    set(fig, 'units','Normalized','Position',windowsize_tests);
        subplot(1,2,1);
            histogram(test_distr_wald_sorted, 'BinWidth', 2, 'Normalization', 'probability');
            hold on;
            histogram(lambda_hat_full, 'BinWidth', 2, 'Normalization', 'probability');
            hold off;
        subplot(1,3,2);
            histogram(test_distr_max_sorted, 'BinWidth', .25, 'Normalization', 'probability');
            hold on;
            histogram(lambda_hat_pars, 'BinWidth', .25, 'Normalization', 'probability');
            hold off;
        suptitle('Non-Linear Model with Weak Identification - Lambdas');
%         temp_title = sprintf('$$ \\hat{\\rho}_{n} $$, %s', innovation_type_string);
%         title(temp_title, 'Interpreter', 'latex', 'FontSize', 22);
%         axis([0, 16, -.4, .4]);
        outputname=sprintf('./%s/%s_J%d_n%d_kl%d', outputdir, 'nonlinear_tests_this_paper', J, n, k_lambda);
%         saveas(fig,outputname,'png');
        
        
        
        
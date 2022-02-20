function time_remaining(j, J)
    if j == 1
        tic;
    end
    if j > 1
        if mod(j,10)==0
            time_remaining_cur = (J - j + 1) * toc / (j-1);
            time_remaining_min = floor(time_remaining_cur/60);
            time_remaining_sec = floor(time_remaining_cur-time_remaining_min*60);
            fprintf('------- \n It: %d/%d \n Est Time = %d min, %d sec \n \n', j, J, time_remaining_min, time_remaining_sec);
        end
    end
end
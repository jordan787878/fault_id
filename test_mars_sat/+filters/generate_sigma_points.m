function [sigma_points, weights, weights_cov] = generate_sigma_points(mu_old, cov_old)
    nx = size(mu_old, 1);
    lambda = 5.0;
    try
        L = chol((nx+lambda)*cov_old, 'lower');
    catch ME
        disp('Matrix is not positive definite.');
        disp(cov_old)
    end
    
    sigma_points = [mu_old];
    weights = [lambda/(nx+lambda)];
    for i = 1:nx
        sigma_points(:,end+1) = mu_old + L(:,i);
        weights(end+1) = 1/(2*(nx+lambda));
    end
    for i = 1:nx
        sigma_points(:,end+1) = mu_old - L(:,i);
        weights(end+1) = 1/(2*(nx+lambda));
    end

    % [tmp]
    weights_cov = weights;
end
function [mu_pri, cov_pri, mu_pos, cov_pos, y_pri, cov_S] = ukf(Y, mu_old, cov_old, params)

mu_pri = 0.0*mu_old; 
L = params.L; Q = params.Q;
cov_pri = L*Q*L';
y_pri = 0.0*Y;

% generate sigma points
[sigma_points, weights] = generate_sigma_points(mu_old, cov_old);

% time update
sigma_points_new = [];
for i = 1:length(weights)
    updated_sigma_point= filters.time_update_state(sigma_points(:,i), params);
    mu_pri = mu_pri + weights(i) * updated_sigma_point;
    sigma_points_new(:,end+1) = updated_sigma_point;
end
for i = 1:length(weights)
    cov_pri = cov_pri + weights(i)*(sigma_points_new(:,i) - mu_pri)*(sigma_points_new(:,i) - mu_pri)';
    y_pri = y_pri + weights(i) * params.H * sigma_points_new(:,i);
end

% measurement update
nx = size(mu_old, 1);
ny = size(Y, 1);
cov_xy = zeros(nx, ny);
cov_yy = zeros(ny, ny);
for i = 1:length(weights)
    cov_xy = cov_xy + weights(i)*(sigma_points_new(:,i) - mu_pri)*(params.H * sigma_points_new(:,i) - y_pri)';
    cov_yy = cov_yy + weights(i)*(params.H * sigma_points_new(:,i) - y_pri)*(params.H * sigma_points_new(:,i) - y_pri)';
end
S = params.R + cov_yy;
cov_S = S;
K = cov_xy*inv(S);
mu_pos = mu_pri + K*(Y - y_pri);
cov_pos = cov_pri - K*S*K';

% Compute eigenvalues
eigenvalues = eig(cov_pri);
% Check if any eigenvalue is complex
if ~isreal(eigenvalues)
    disp('cov_pri has complex eigenvalues and is not PSD.');
    disp(cov_pri)
    disp(sigma_points_new)
% Check if any eigenvalue is significantly negative
elseif any(eigenvalues < 0)
    disp('cov_pri is not PSD due to negative eigenvalues.');
    disp(cov_pri)
    disp(sigma_points)
end

end

function [sigma_points, weights] = generate_sigma_points(mu_old, cov_old)
nx = size(mu_old, 1);
% kappa = 0; alpha = 10^-4;
% lambda = alpha^2*(nx+kappa) - nx;
lambda = 1.0;

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
end
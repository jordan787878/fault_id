function u_active = active_control(sim, t1, t2, mu_H_pos, cov_H_pos, Hypothesis)
% Description
%   Note that the separation metric is weighted by prior beliefs
%   From the linear approximation, the optimal solutions are the vertices of the control set.
%   This makes it efficient to determine optimal control
%   We also try the mutual information metric, yet it is expensive to
%   compute due to the integral over all possible measurements.

nu = 3;

% Define safe boundaries for each component (3x1 vectors)
lb = -1e-3*ones(1,nu);  % Lower bounds for each row
ub =  1e-3*ones(1,nu);   % Upper bounds for each row
% lb = -5e-4*ones(1,nu);  % Lower bounds for each row
% ub =  5e-4*ones(1,nu);   % Upper bounds for each row

safe_lower = lb;
safe_upper = ub;
% Generate all combinations (vertex coordinates) using ndgrid:
[X, Y, Z] = ndgrid([safe_lower(1), safe_upper(1)], ...
                   [safe_lower(2), safe_upper(2)], ...
                   [safe_lower(3), safe_upper(3)]);

% Combine them into a 3x8 matrix where each column is a vertex vector:
U_safe = [X(:)'; Y(:)'; Z(:)'];
U_safe(:,end+1) = [0;0;0];

% Number of random samples to generate
% N = 100;  % Adjust as needed
% % Preallocate U_safe (nu x N)
% U_safe = zeros(nu, N);
% % Generate U_safe: For each component i, sample uniformly between safe_lower(i) and safe_upper(i)
% for i = 1:nu
%     U_safe(i, :) = safe_lower(i) + (safe_upper(i) - safe_lower(i)) * rand(1, N);
% end

separate_metrics = [];
I_metrics = [];
for i = 1:size(U_safe,2)
    separate_metrics(end+1) = get_separate_metric(sim, t1, t2, mu_H_pos, cov_H_pos, U_safe(:,i), Hypothesis);

    % I_metrics(end+1) = get_mutual_information(sim, t1, t2, mu_H_pos, cov_H_pos, U_safe(:,i));
end

% disp(I_metrics)
% [~, index_I] = max(I_metrics);
% disp(num2str(index) + " , " +num2str(index_I))
% u_active = U_safe(index_I);

% disp(separate_metrics)
[~, index] = max(separate_metrics);
u_active = U_safe(index);
end


function metric = get_separate_metric(sim, t1, t2, mu_H_pos, cov_H_pos, u, Hypothesis)
metric = 0;
min_dist = inf;
geo_mean = 1;
n_Hypo = size(mu_H_pos,2)-1;
for m_i = 0:n_Hypo
    for m_j = 0:n_Hypo
        if(m_i == m_j)
            metric = metric + 0;
        else
            hypo_i = Hypothesis(m_i+1);
            hypo_j = Hypothesis(m_j+1);
            [mu_y_i, cov_y_i]       = wrapper_ukf(sim, t1, t2, mu_H_pos(:,m_i+1), cov_H_pos(:,:,m_i+1), u, m_i);
            [mu_y_j, cov_y_j] = wrapper_ukf(sim, t1, t2, mu_H_pos(:,m_j+1), cov_H_pos(:,:,m_j+1), u, m_j);
            distance = get_Mah_distance(mu_y_i, mu_y_j, cov_y_j)*sqrt(hypo_i)*sqrt(hypo_j);
            metric = metric + distance;
            min_dist = min(min_dist, distance);
            geo_mean = geo_mean * distance;
        end
    end
end
% metric = metric * min_dist;
metric = geo_mean^(1/(n_Hypo^2));
% metric = min_dist;
end


function [y_pri, cov_S] = wrapper_ukf(sim, t1, t2, mu_old, cov_old, u, mode)

params = dynamics.get_mode_dyn_params(u, sim, mode);
mu_pri = 0.0*mu_old; 
L = params.L; Q = params.Q; R = params.R;
ny = size(R, 1);
cov_pri = L*Q*L';
y_pri = zeros(ny, 1);
% generate sigma points
[sigma_points, weights] = filters.generate_sigma_points(mu_old, cov_old);
% time update
sigma_points_new = [];
for i = 1:length(weights)
    updated_sigma_point = filters.time_update_state(sim, sigma_points(:,i), u, t1, t2, mode);
    mu_pri = mu_pri + weights(i) * updated_sigma_point;
    sigma_points_new(:,end+1) = updated_sigma_point;
end
for i = 1:length(weights)
    cov_pri = cov_pri + weights(i)*(sigma_points_new(:,i) - mu_pri)*(sigma_points_new(:,i) - mu_pri)';
    y_pri = y_pri + weights(i) * params.H * sigma_points_new(:,i);
end
cov_yy = zeros(ny, ny);
for i = 1:length(weights)
    cov_yy = cov_yy + weights(i)*(params.H * sigma_points_new(:,i) - y_pri)*(params.H * sigma_points_new(:,i) - y_pri)';
end
cov_S = R + cov_yy; 

end


function d = get_Mah_distance(x, mu, cov)
d = sqrt((x-mu)'*cov^(-1)*(x-mu));
end


% function I = get_mutual_information(sim, t1, t2, mu_H_pos, cov_H_pos, u)
% n_Hypo = size(mu_H_pos,2);
% mu_y = zeros(6, n_Hypo);
% cov_y = zeros(6,6,n_Hypo);
% 
% for m = 0:n_Hypo-1
%     [mu_y_i, cov_y_i] = wrapper_ukf(sim, t1, t2, mu_H_pos(:,m+1), cov_H_pos(:,:,m+1), u, m);
%     mu_y(:,m+1) = mu_y_i;
%     cov_y(:,:,m+1) = cov_y_i;
% end
% 
% Y_1 = sum(mu_y, 2);
% 
% likelihoods = [];
% for m = 0:n_Hypo-1
%     % like = evaluate_likelihood(Y_1, mu_y(:, m+1), cov_y(:, :, m+1));
%     % likelihoods(end+1) = like;
% 
%     S_used = cov_y(:,:,m+1);
%     inno_used = Y_1 - mu_y(:, m+1);
%     L_lower = chol(S_used, 'lower');
%     % Compute log(det(S_used)) as 2 * sum(log(diag(L)))
%     logDet_S = 2 * sum(log(diag(L_lower)));
%     % Solve L * x = inno_used to compute the quadratic form
%     x = L_lower \ inno_used;
%     quadForm = sum(x.^2);
%     % Compute the log likelihood
%     log_like = -0.5 * ( length(inno_used) * log(2*pi) + logDet_S + quadForm )
%     like = exp(log_like);
%     likelihoods(end+1) = like;
% end
% disp("likes " + num2str(likelihoods))
% sum_likelihoods = sum(likelihoods);
% 
% I = 0;
% for m = 0:n_Hypo-1
%     I = I + likelihoods(m+1)*log(likelihoods(m+1)/sum_likelihoods);
% end
% end
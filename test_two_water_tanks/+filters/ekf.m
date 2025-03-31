function [mu_pri, cov_pri, mu_pos, cov_pos, y_pri, cov_S] = ekf(Y, mu_old, cov_old, params)
[mu_pri, cov_pri] = filters.ekf_time_update(mu_old, cov_old, params);
y_pri = params.H * mu_pri;
cov_S = params.H*cov_pri*params.H' + params.R;
H = params.H;
R = params.R;
K = cov_pri*H'*(H*cov_pri*H' + R)^(-1);
mu_pos = mu_pri + K*(Y - H*mu_pri);
cov_pos = (eye(length(R))-K*H)*cov_pri*(eye(length(R))-K*H)' + K*R*K';
end
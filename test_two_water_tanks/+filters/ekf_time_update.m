function [mu, cov] = ekf_time_update(mu_old, cov_old, params)
    mu = filters.time_update_state(mu_old, params);
    F = predict_Jacobian_F(mu_old, params);
    L = params.L; Q = params.Q;
    cov = F*cov_old*F' + L*Q*L';
end

function F = predict_Jacobian_F(mu_old, params)
% u = params.u;
C1L = params.C1L;
C12 = params.C12;
C2L = params.C2L;
dt  = params.dt;
A   = params.A;
x1 = mu_old(1);
x2 = mu_old(2);
F = zeros(2,2);
F(1,1) = 1 - (dt*(C1L/(2*x1^(1/2)) + C12/(2*(x1 - x2)^(1/2))))/A;
F(1,2) = (C12*dt)/(2*A*(x1 - x2)^(1/2));
F(2,1) = (C12*dt)/(2*A*(x1 - x2)^(1/2));
F(2,2) = 1 - (dt*(C2L/(2*x2^(1/2)) + C12/(2*(x1 - x2)^(1/2))))/A;
end
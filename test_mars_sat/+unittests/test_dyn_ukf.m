function [sim, fid] = test_dyn_ukf()
rng(0)

sim.constants.cov_i = 1e-3*eye(6);
sim.constants.L_ext = [0;0;0];
sim.constants.I_inertia = diag([10, 5, 7.5]);
sim.constants.L = [zeros(3,3); eye(3)];
sim.constants.Q = 1e-9*eye(3);
sim.constants.H = eye(6);
sim.constants.R = 2e-5*eye(6);
sim.constants.controller_K = (1/6)^2/5;
sim.constants.controller_P = 1/6;
sim.config.true_mode = 0;

start_time = 0;
inc_time = 1;
end_time = 7000;

% initialization
nx = size(sim.constants.L, 1);
ny = size(sim.constants.H, 1);
nw = size(sim.constants.Q, 1);
nv = size(sim.constants.R, 1);
t_hist = start_time:inc_time:end_time;
x_hist = [satsym.set_x_init(start_time)];
y_hist = [zeros(ny,1)];
[ctrl_torque, ctrl_mode, r_LMO, r_GMO, BN, w_BR_Body] = satsym.get_control_and_mode(t_hist(end), x_hist(:,end), 0, sim);
rlmo_hist = [r_LMO]; rgmo_hist = [r_GMO]; bn_hist = [BN]; ctrlmode_hist = [ctrl_mode];
u_hist = [ctrl_torque];

% 4 Hypothesis
n_Hypo = 4;
Hypothesis = ones(1, n_Hypo)*(1/(n_Hypo)); 
x_i_pos = x_hist(:,end);
mu_H_pos = repmat(x_i_pos, [1, n_Hypo]);
cov_H_pos = repmat(sim.constants.cov_i, [1, 1, n_Hypo]);
mu_H_pri = repmat(x_i_pos, [1, n_Hypo]);
cov_H_pri = repmat(sim.constants.cov_i, [1, 1, n_Hypo]);
H_hist = Hypothesis;
mu_H_pos_hist = mu_H_pos;
cov_H_pos_hist = cov_H_pos;
mu_H_pri_hist = mu_H_pri;
cov_H_pri_hist = cov_H_pri;
inno_H_pri = repmat(zeros(ny,1), [1, n_Hypo]);
inno_H_pri_hist = inno_H_pri;
S_H_pri = repmat(sim.constants.H*sim.constants.cov_i*sim.constants.H', [1, 1, n_Hypo]);
S_H_pri_hist = S_H_pri;

for k = 1:length(t_hist)-1
    t1 = t_hist(k); t2 = t_hist(k+1);
    % propagate from t1 to t2
    w = mvnrnd(zeros(nw, 1), sim.constants.Q, 1)';
    x_hist(:,end+1) = satsym.update_state(sim, x_hist(:,end), ctrl_torque, t1, t2, w, sim.config.true_mode);
    % measurement at t2
    v = mvnrnd(zeros(nv,1), sim.constants.R, 1)';
    y_hist(:,end+1) = sim.constants.H*x_hist(:,end) + v;

    % Filtering
    Y = y_hist(:,end);
    u_old = ctrl_torque;
    for m = 0:n_Hypo-1
        mu_old = mu_H_pos(:,m+1); cov_old = cov_H_pos(:,:,m+1);
    
        [mu_pri, cov_pri, mu_pos, cov_pos, y_pri, cov_S] = filters.ukf(sim, t1, t2, Y, mu_old, cov_old, u_old, m);
    
        mu_H_pos(:,m+1) = mu_pos; cov_H_pos(:, :, m+1) = cov_pos;
        inno_H_pri(:, m+1) = Y - y_pri; S_H_pri(:,:,m+1) = cov_S;
    end

    % determine control at t2
    [ctrl_torque, ctrl_mode, r_LMO, r_GMO, BN, w_BR_Body] = satsym.get_control_and_mode(t2, x_hist(:,end), 0, sim);

    % store data
    u_hist(:, end+1) = ctrl_torque;
    rlmo_hist(:, end+1) = r_LMO;
    rgmo_hist(:, end+1) = r_GMO;
    ctrlmode_hist(end+1) = ctrl_mode;
    bn_hist(:,:,end+1) = BN;
    mu_H_pos_hist = vertcat(mu_H_pos_hist, mu_H_pos);
    cov_H_pos_hist = vertcat(cov_H_pos_hist, cov_H_pos);
    inno_H_pri_hist = vertcat(inno_H_pri_hist, inno_H_pri);
    S_H_pri_hist = vertcat(S_H_pri_hist, S_H_pri);
end

% prepare output data
fid.outputs.t_hist = t_hist;
fid.outputs.x_hist = x_hist;
fid.outputs.y_hist = y_hist;
fid.outputs.u_hist = u_hist;
fid.outputs.rlmo_hist = rlmo_hist;
fid.outputs.rgmo_hist = rgmo_hist;
fid.outputs.bn_hist = bn_hist;
fid.outputs.ctrlmode_hist = ctrlmode_hist;
fid.outputs.mu_H_pos_hist = mu_H_pos_hist;
fid.outputs.cov_H_pos_hist = cov_H_pos_hist;
fid.outputs.inno_H_pri_hist = inno_H_pri_hist;
fid.outputs.S_H_pri_hist = S_H_pri_hist;

end
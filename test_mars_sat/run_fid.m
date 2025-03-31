function [sim, fid] = run_fid(sim, fid, true_mode)

% Description
%   This scripts uses satsym to demonstrate active fault identification
%   The simulation will run until the end time; the fault id will run until
%   termination.
%   The sim stores the constants and configurations of the simulation
%   The fid stores the configurations of the FID manager, and the results in fid.outputs

k_stop = 25 + fid.config.moving_window*2.0;

start_time = sim.config.start_time;
inc_time = sim.config.inc_time;
end_time = sim.config.end_time;

% if doing multi trials, direct end simulation once fid stops running
if(sim.config.N_trials > 1)
    end_time = start_time + k_stop + 1;
end

% initialization
nx = size(sim.constants.L, 1);
ny = size(sim.constants.H, 1);
nw = size(sim.constants.Q, 1);
nv = size(sim.constants.R, 1);
t_span = start_time:inc_time:end_time;
t_hist = [t_span(1)];

x_i_nominal = satsym.set_x_init(start_time);
x_i = x_i_nominal + mvnrnd(zeros(nx, 1), sim.constants.cov_i, 1)';
x_hist = [x_i];

y_hist = [zeros(ny,1)];
[ctrl_torque, ctrl_mode, r_LMO, r_GMO, BN, w_BR_Body] = satsym.get_control_and_mode(t_span(1), x_hist(:,end), 0, sim);
rlmo_hist = [r_LMO]; rgmo_hist = [r_GMO]; bn_hist = [BN]; ctrlmode_hist = [ctrl_mode];
u_hist = [ctrl_torque];
running_fid_flag = true;

fid.outputs.true_mode = true_mode;
fid.outputs.id_mode = NaN;
fid.outputs.fail   = 0; df = ny;
chi_crit = chi2inv(1-fid.config.alpha, fid.config.moving_window*ny)/fid.config.moving_window;
chi_crit_up = chi2inv(1-fid.config.alpha/2, df*fid.config.moving_window)/fid.config.moving_window;
chi_crit_lo = chi2inv(fid.config.alpha/2, df*fid.config.moving_window)/fid.config.moving_window;

% 4 Hypothesis
n_Hypo = fid.config.N_hypo;
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

for k = 1:length(t_span)-1
    t1 = t_span(k); t2 = t_span(k+1);
    u_old = u_hist(:,end);
    % propagate from t1 to t2
    w = mvnrnd(zeros(nw, 1), sim.constants.Q, 1)';
    x_hist(:,end+1) = satsym.update_state(sim, x_hist(:,end), u_old, t1, t2, w, true_mode);
    % measurement at t2
    v = mvnrnd(zeros(nv,1), sim.constants.R, 1)';
    y_hist(:,end+1) = sim.constants.H*x_hist(:,end) + v;

    % Filtering
    if(running_fid_flag)
    Y = y_hist(:,end);
    for m = 0:n_Hypo-1
        mu_old = mu_H_pos(:,m+1); cov_old = cov_H_pos(:,:,m+1);
    
        [mu_pri, cov_pri, mu_pos, cov_pos, y_pri, cov_S] = filters.ukf(sim, t1, t2, Y, mu_old, cov_old, u_old, m);
    
        mu_H_pos(:,m+1) = mu_pos; cov_H_pos(:, :, m+1) = cov_pos;
        inno_H_pri(:, m+1) = Y - y_pri; S_H_pri(:,:,m+1) = cov_S;
    end
    % store filtering data
    mu_H_pos_hist = vertcat(mu_H_pos_hist, mu_H_pos);
    cov_H_pos_hist = vertcat(cov_H_pos_hist, cov_H_pos);
    inno_H_pri_hist = vertcat(inno_H_pri_hist, inno_H_pri);
    S_H_pri_hist = vertcat(S_H_pri_hist, S_H_pri);

    % Hypothesis update
    if(k >= fid.config.moving_window)
        % likelihood of each hypothesis
        logL = zeros(1,n_Hypo);
        for m = 0:n_Hypo-1
            inno_used = inno_H_pri_hist(end-(ny*fid.config.moving_window-1):end, m+1);
            tmp = S_H_pri_hist(end-(ny*fid.config.moving_window-1):end, : , m+1);
            S_used = zeros(ny*fid.config.moving_window, ny*fid.config.moving_window);
            for j = 1:fid.config.moving_window
                S_used(ny*(j-1)+1:ny*(j-1)+ny, ny*(j-1)+1:ny*(j-1)+ny) = tmp(ny*(j-1)+1:ny*(j-1)+ny, :);
            end
            chi = (inno_used)'*inv(S_used)*(inno_used)/fid.config.moving_window;

            inno_whiten = reshape(inno_used, ny, fid.config.moving_window);
            for jj = 0:fid.config.moving_window-1
                inno_jj = inno_whiten(:,jj+1);
                S_jj   = S_used(ny*jj+1:ny*jj+ny, ny*jj+1:ny*jj+ny);
                inno_whiten_jj = inv(chol(S_jj, 'lower'))*inno_jj;
                inno_whiten(:,jj+1) = inno_whiten_jj;
            end
            chi_sr = sum(inno_whiten(:).^2)/fid.config.moving_window;
            % disp(chi - chi_sr)
            inno_used_whiten = reshape(inno_whiten, ny*fid.config.moving_window, 1);
            S_used_whiten = eye(ny*fid.config.moving_window);
            
            % if(chi_sr < chi_crit)
            if(chi_sr <= chi_crit_up && chi_sr >= chi_crit_lo)
            
                % likelihood_y = filters.evaluate_likelihood(inno_used_whiten, 0.0*inno_used_whiten, S_used_whiten);
                % p_post = likelihood_y * Hypothesis(m+1);
                % Hypothesis(m+1) = p_post;
                logP = -0.5*log_determinant(S_used) -0.5*chi_sr*fid.config.moving_window + log(Hypothesis(m+1));
                Hypothesis(m+1) = logP;
            else
                % Hypothesis(m+1) = 0;
                Hypothesis(m+1) = -Inf;
            end
        end
        % if(sum(Hypothesis) == 0)
        %     % all hypothesis are rejected - do renormalization
        %     Hypothesis = ones(1, n_Hypo)*(1/(n_Hypo));
        % else
        %     Hypothesis = Hypothesis/(sum(Hypothesis));
        % end
        maxlogP = max(Hypothesis);
        if(maxlogP == -Inf)
            Hypothesis = ones(1, n_Hypo)*(1/(n_Hypo));
        else
            indices_not_minusInf = find(Hypothesis ~= -Inf);
            logP_not_minusInf = Hypothesis(indices_not_minusInf);
            P_prime = exp(logP_not_minusInf - maxlogP);
            sum_P_prime = sum(P_prime);
            P_subset = P_prime/sum_P_prime;
            Hypothesis(indices_not_minusInf) = P_subset;
            Hypothesis(Hypothesis == -Inf) = 0;
            % disp(Hypothesis)
            % disp("==")
        end
    end
    % store hypothesis history % disp("Hypothesis at t: " +num2str(t2) + ", " + num2str(Hypothesis))
    H_hist = vertcat(H_hist, Hypothesis);
    end

    % Control at t2
    [u_apply, ctrl_mode, r_LMO, r_GMO, BN, ~] = satsym.get_control_and_mode(t2, x_hist(:,end), 0, sim);
    if(running_fid_flag && fid.config.action_flag)
        u_apply = active_control(sim, t2, t2+sim.config.inc_time, mu_H_pos, cov_H_pos, Hypothesis);
    end

    % store control data
    t_hist(end+1) = t2;
    u_hist(:, end+1) = u_apply;
    rlmo_hist(:, end+1) = r_LMO;
    rgmo_hist(:, end+1) = r_GMO;
    ctrlmode_hist(end+1) = ctrl_mode;
    bn_hist(:,:,end+1) = BN;

    % Termination condition
    if(running_fid_flag)
        [terminate_flag, k_end, fail, id_mode] = event_termination(true_mode, k, k_stop, fid, Hypothesis);
        if(terminate_flag)
            running_fid_flag = false;
            fid.outputs.k_end = k_end+1; fid.outputs.fail = fail; 
            fid.outputs.id_mode = id_mode;
        end
    end
end

fid.outputs.t_hist = t_hist;
fid.outputs.x_hist = x_hist;
fid.outputs.y_hist = y_hist;
fid.outputs.u_hist = u_hist;
fid.outputs.H_hist = H_hist;
fid.outputs.rlmo_hist = rlmo_hist;
fid.outputs.rgmo_hist = rgmo_hist;
fid.outputs.bn_hist = bn_hist;
fid.outputs.ctrlmode_hist = ctrlmode_hist;
fid.outputs.mu_H_pos_hist = mu_H_pos_hist;
fid.outputs.cov_H_pos_hist = cov_H_pos_hist;
fid.outputs.inno_H_pri_hist = inno_H_pri_hist;
fid.outputs.S_H_pri_hist = S_H_pri_hist;

end


function [terminate_flag, k_end, fail, id_mode] = event_termination(true_mode, k, k_stop, fid, Hypothesis)
terminate_flag = false; k_end = NaN; fail  = NaN; id_mode = NaN;
if(k == k_stop)
    if(true_mode == -1)
        if(max(Hypothesis) < fid.config.crit)
            k_end = k;  % count
            fail = 0;
            % fid.outputs.k_end = -1; % does not count this for unknown faults
            % disp("Success, unknown fault identified")
            id_mode = -1;
            terminate_flag = true;
            return;
        end
    else
        if(max(Hypothesis) < fid.config.crit)
            k_end = k;
            fail = 1;
            % disp("...Fail, true:"+num2str(true_mode)+", not detected over kf: " +num2str(k_stop))
            terminate_flag = true;
            return;
        end
    end
end
% 2. false hypothesis is identified
false_Hypothesis = Hypothesis([1:true_mode+1-1, true_mode+1+1:end]);
if(max(false_Hypothesis) > fid.config.crit)
    [~, detected] = max(Hypothesis);
    % disp("...Fail, true:"+num2str(true_mode)+", detected:"+num2str(detected-1))
    id_mode = detected-1;
    k_end = k;
    fail = 1;
    terminate_flag = true;
    return;
end
% 3. true hypothesis is identified
if(true_mode+1 > 0 && Hypothesis(true_mode+1) > fid.config.crit)
    % disp("Success, detect fault "+num2str(true_mode))
    [~,index] = max(Hypothesis);
    id_mode = index-1;
    k_end = k;
    fail = 0;
    terminate_flag = true;
    return;
end
end

function logDet = log_determinant(S)
    [L, U, P] = lu(S);
    diagU = diag(U);
    signU = prod(sign(diagU));
    logDet = sum(log(abs(diagU)));
    if signU < 0
        logDet = logDet + log(-1); % Account for negative determinant
    end
end
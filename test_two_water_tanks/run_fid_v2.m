function [sim, fid] = run_fid_v2(sim, fid, true_mode)

% Description
%   add belief weights to separation metric in active control
% disp(true_mode)

constants = sim.constants;

% sim initialization
x_i = mvnrnd(constants.mu_i, constants.cov_i, 1)';
while(x_i(1) - 0.5 < x_i(2) + 0.5 || x_i(2) < fid.config.x2_min + 1.0 || x_i(1) + 1.0 > fid.config.x1_max)
    x_i = mvnrnd(constants.mu_i, constants.cov_i, 1)';
end
x_hist = [x_i];

% fid initialization
fid.outputs.true_mode = true_mode;
fid.outputs.id_mode = NaN;
fid.outputs.k_end = inf;
fid.outputs.fail   = 0;
ny = size(sim.constants.R, 2); df = ny;
chi_crit = chi2inv(1-fid.config.alpha, fid.config.moving_window*df)/fid.config.moving_window;
chi_crit_up = chi2inv(1-fid.config.alpha/2, df*fid.config.moving_window)/fid.config.moving_window;
chi_crit_lo = chi2inv(fid.config.alpha/2, df*fid.config.moving_window)/fid.config.moving_window;
k_stop = 20 + fid.config.moving_window*2.0;
fid.outputs.u_hist = [fid.config.u_default];
fid.outputs.y_hist = [zeros(2,1)];
fid.outputs.M_hist = [];

% 4 Hypothesis
n_Hypo = fid.config.N_hypo;
Hypothesis = ones(1, n_Hypo)*(1/(n_Hypo));
x_i_pos = constants.mu_i; 
mu_H_pos = repmat(x_i_pos, [1, n_Hypo]);
cov_H_pos = repmat(constants.cov_i, [1, 1, n_Hypo]);
mu_H_pri = repmat(x_i_pos, [1, n_Hypo]);
cov_H_pri = repmat(constants.cov_i, [1, 1, n_Hypo]);
H_hist = Hypothesis;
mu_H_pos_hist = mu_H_pos;
cov_H_pos_hist = cov_H_pos;
mu_H_pri_hist = mu_H_pri;
cov_H_pri_hist = cov_H_pri;
inno_H_pri = repmat(zeros(2,1), [1, n_Hypo]);
inno_H_pri_hist = inno_H_pri;
S_H_pri = repmat(constants.H*constants.cov_i*constants.H', [1, 1, n_Hypo]);
S_H_pri_hist = S_H_pri;
M_U = [];

for k = 1:k_stop
    % Get current true state and measurement
    if(sim.config.use_process_noise)
        w = mvnrnd([0;0], constants.Q, 1)';
    else
        w = [0;0];
    end
    x_hist(:,end+1) = dynamics.update_state(x_hist(:,end), fid.outputs.u_hist(:,end), w, constants, true_mode);
    if(sim.config.use_measure_noise)
        v = mvnrnd([0;0], constants.R, 1)';
        fid.outputs.y_hist(:,end+1) = x_hist(:,end) + v;
    else
        fid.outputs.y_hist(:,end+1) = x_hist(:,end);
    end

    % Filtering
    Y = fid.outputs.y_hist(:, end);
    ny = size(Y,1);
    u_old = fid.outputs.u_hist(:,end);
    for m = 0:n_Hypo-1
        mu_old = mu_H_pos(:,m+1); cov_old = cov_H_pos(:,:,m+1);
        params = dynamics.get_mode_dyn_params(u_old, constants, m);

        if(fid.config.filter_options == "ukf")
            [mu_pri, cov_pri, mu_pos, cov_pos, y_pri, cov_S] = filters.ukf(Y, mu_old, cov_old, params); % ukf
        else
            [mu_pri, cov_pri, mu_pos, cov_pos, y_pri, cov_S] = filters.ekf(Y, mu_old, cov_old, params); % efk
        end

        mu_H_pri(:,m+1) = mu_pri; cov_H_pri(:,:,m+1) = cov_pri;
        mu_H_pos(:,m+1) = mu_pos; cov_H_pos(:, :, m+1) = cov_pos;
        inno_H_pri(:, m+1) = Y - y_pri; S_H_pri(:,:,m+1) = cov_S;
    end
    % store pri and pos uncertainty estimates
    mu_H_pri_hist = vertcat(mu_H_pri_hist, mu_H_pri);
    cov_H_pri_hist = vertcat(cov_H_pri_hist, cov_H_pri);
    mu_H_pos_hist = vertcat(mu_H_pos_hist, mu_H_pos);
    cov_H_pos_hist = vertcat(cov_H_pos_hist, cov_H_pos);
    inno_H_pri_hist = vertcat(inno_H_pri_hist, inno_H_pri);
    S_H_pri_hist = vertcat(S_H_pri_hist, S_H_pri);

    % Hypothesis update
    if(k >= fid.config.moving_window)
        for m = 0:n_Hypo-1
            inno_used = inno_H_pri_hist(end-(ny*fid.config.moving_window-1):end, m+1);
            tmp = S_H_pri_hist(end-(ny*fid.config.moving_window-1):end, : , m+1);
            S_used = zeros(ny*fid.config.moving_window, ny*fid.config.moving_window);
            for j = 1:fid.config.moving_window
                S_used(ny*j-1:ny*j, ny*j-1:ny*j) = tmp(ny*j-1:ny*j, :);
            end
            if(fid.config.likeli_options == "1")
                chi = (inno_used)'*inv(S_used)*(inno_used)/fid.config.moving_window;
            else
                inno_whiten = reshape(inno_used, ny, fid.config.moving_window);
                for jj = 0:fid.config.moving_window-1
                    inno_jj = inno_whiten(:,jj+1);
                    S_jj   = S_used(ny*jj+1:ny*jj+ny, ny*jj+1:ny*jj+ny);
                    L_jj = inv(chol(S_jj, 'lower'));
                    check_value = L_jj'*L_jj - inv(S_jj);
                    if(abs(check_value) > 1e-10)
                        disp("check "+num2str(check_value))
                    end
                    inno_whiten_jj = L_jj*inno_jj;
                    inno_whiten(:,jj+1) = inno_whiten_jj;
                end
                chi = sum(inno_whiten(:).^2)/fid.config.moving_window;
                % inno_used_whiten = reshape(inno_whiten, ny*fid.config.moving_window, 1);
                % S_used_whiten = eye(ny*fid.config.moving_window);
            end

            % check
            % disp(num2str(m)+" "+num2str(chi)+" "+num2str(chi_crit)+" "+num2str(chi_crit_lo)+" "+num2str(chi_crit_up))
            if(fid.config.reject_flag)
                % (option 1) with hypothesis rejecting process
                % if(chi < chi_crit) % (option 1.1 - one-tailed test)
                if(chi <= chi_crit_up && chi >= chi_crit_lo) % (option 1.2 - two-tailed test)
                    % if(fid.config.likeli_options == "1")
                    %     likelihood_y = filters.evaluate_likelihood(inno_used, 0.0*inno_used, S_used);
                    % else
                    %     likelihood_y = filters.evaluate_likelihood(inno_used_whiten, 0.0*inno_used_whiten, S_used_whiten);
                    % end
                    % p_post = likelihood_y * Hypothesis(m+1);
                    % Hypothesis(m+1) = p_post;
                    logP = -0.5*log_determinant(S_used) -0.5*chi*fid.config.moving_window + log(Hypothesis(m+1));
                    Hypothesis(m+1) = logP;
                else
                    % Hypothesis(m+1) = 0;
                    Hypothesis(m+1) = -Inf;
                end
            else
                logP = -0.5*log_determinant(S_used) -0.5*chi*fid.config.moving_window + log(Hypothesis(m+1));
                Hypothesis(m+1) = logP;
            end

        end
        % (option 1 - do renormalization)
        if(fid.config.renorm_flag)
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
            end
            % % normalize posterior probability
            % if(sum(Hypothesis) == 0.0)
            %     % all hypothesis are rejected - do renormalization 
            %     Hypothesis = ones(1, n_Hypo)*(1/(n_Hypo));
            % else
            %     Hypothesis = Hypothesis/(sum(Hypothesis));
            % end
        % (option 2 - no renorm)
        else
            maxlogP = max(Hypothesis);
            if(maxlogP == -Inf)
                Hypothesis = zeros(1, n_Hypo);
            else
                indices_not_minusInf = find(Hypothesis ~= -Inf);
                logP_not_minusInf = Hypothesis(indices_not_minusInf);
                P_prime = exp(logP_not_minusInf - maxlogP);
                sum_P_prime = sum(P_prime);
                P_subset = P_prime/sum_P_prime;
                Hypothesis(indices_not_minusInf) = P_subset;
                Hypothesis(Hypothesis == -Inf) = 0;
            end
        end
    end
    % store hypothesis history
    H_hist = vertcat(H_hist, Hypothesis);

    % Active control
    if(fid.config.action_flag)
        u_apply = active_control(fid, mu_H_pos, cov_H_pos, constants, Hypothesis);
    else
        u_apply = fid.config.u_default;
    end
    % store data
    fid.outputs.u_hist(:,end+1) = u_apply;
    
    % Termination condition
    [terminate_flag, k_end, fail, id_mode] = event_termination(true_mode, k, k_stop, fid, Hypothesis, x_hist(:,end));
    if(terminate_flag)
        fid.outputs.k_end = k_end; fid.outputs.fail = fail; 
        fid.outputs.mu_H_pri_hist = mu_H_pri_hist;
        fid.outputs.cov_H_pri_hist = cov_H_pri_hist;
        fid.outputs.mu_H_pos_hist = mu_H_pos_hist;
        fid.outputs.cov_H_pos_hist = cov_H_pos_hist;
        fid.outputs.inno_H_pri_hist = inno_H_pri_hist;
        fid.outputs.S_H_pri_hist = S_H_pri_hist;
        fid.outputs.H_hist = H_hist;
        fid.outputs.x_hist = x_hist;
        fid.outputs.id_mode = id_mode;
        break;
    end
end

end


% [updated with return;]
function [terminate_flag, k_end, fail, id_mode] = event_termination(true_mode, k, k_stop, fid, Hypothesis, x)
terminate_flag = false; k_end = NaN; fail  = NaN; id_mode = NaN;

% if(~isreal(x))
%     terminate_flag = true;
%     fail = 1;
%     disp("FAIL, state constraints violate")
%     return;
% end
if(sum(Hypothesis) == 0.0 && true_mode > -1)
    % disp("FAIL, true mode " + num2str(true_mode) + " all reject ")
    k_end = k;
    fail = 1;
    id_mode = -1;
    terminate_flag = true;
    return;
end

if(k == k_stop)
    if(true_mode == -1)
        if(max(Hypothesis) < fid.config.crit)
            k_end = k;  % count
            fail = 0;
            % fid.outputs.k_end = -1; % does not count this for unknown faults
            % disp("success, unknown fault identified")
            id_mode = -1;
            terminate_flag = true;
            return;
        end
    else
        if(max(Hypothesis) < fid.config.crit)
            k_end = k;
            fail = 1;
            % disp("FAIL, true:"+num2str(true_mode)+", not detected over kf: " +num2str(k_stop))
            terminate_flag = true;
            return;
        end
    end
end
% 2. false hypothesis is identified
false_Hypothesis = Hypothesis([1:true_mode+1-1, true_mode+1+1:end]);
if(max(false_Hypothesis) > fid.config.crit)
    [~, detected] = max(Hypothesis);
    % disp("FAIL, true:"+num2str(true_mode)+", detected:"+num2str(detected-1))
    id_mode = detected-1;
    k_end = k;
    fail = 1;
    terminate_flag = true;
    return;
end
% 3. true hypothesis is identified
if(true_mode+1 > 0 && Hypothesis(true_mode+1) > fid.config.crit)
    % disp("success, detect fault "+num2str(true_mode))
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
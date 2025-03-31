function u_active = active_control(fid, mu_H_pos, cov_H_pos, constants, Hypothesis)

U_safe = get_safecontrol_set(fid, mu_H_pos, cov_H_pos, constants);

if(length(U_safe) == 1)
    u_active = U_safe(end);
else
    separate_metrics = [];
    for i = 1:length(U_safe)
        separate_metrics(end+1) = get_separate_metric(U_safe(i), mu_H_pos, cov_H_pos, constants, Hypothesis);
    end
    if(all(separate_metrics == separate_metrics(1)))
        u_active = fid.config.u_default;
    else
        [~, index] = max(separate_metrics);
        u_active = U_safe(index);
    end
end

end

function U_safe = get_safecontrol_set(fid, mu_H_pos, cov_H_pos, constants)
n_Hypo = length(mu_H_pos)-1;
U_span = 0:0.0001:0.01;
U_candidate = [];
margin = 1.0;
 
for j = 1:length(U_span)
    u_j = U_span(j);
    unsafe_flag = false;
    for m = 0:n_Hypo
        params = dynamics.get_mode_dyn_params(u_j, constants, m);
        [mu_x, cov_x] = filters.ekf_time_update(mu_H_pos(:,m+1), cov_H_pos(:,:,m+1), params);
        if(mu_x(1)+3*cov_x(1,1) > fid.config.x1_max-margin)
           unsafe_flag = true;
           break;
        end
        if(mu_x(2) < fid.config.x2_min + 3*cov_x(2,2) + 1*margin + 0.5) % 1*margin (default)
           unsafe_flag = true;
           break;
        end
        if(mu_x(1)-3*cov_x(1,1)-margin < mu_x(2)+3*cov_x(2,2)+margin)
           unsafe_flag = true;
           break;
        end
    end
    if(unsafe_flag == false)
        U_candidate(:,end+1) = u_j;
    end
end

if(length(U_candidate) < 1)
    U_safe = 0;
elseif(max(U_candidate) < fid.config.u_default)
    U_safe = [0, max(U_candidate)];
else
    U_safe = [0, fid.config.u_default, max(U_candidate)];
end

end

function metric = get_separate_metric(u, mu_H_pos, cov_H_pos, constants, Hypothesis)
metric = 0;
min_dist = inf;
geo_mean = 1;
n_Hypo = length(mu_H_pos)-1;
for m_i = 0:n_Hypo
    for m_j = 0:n_Hypo
        if(m_i == m_j)
            metric = metric + 0;
        else
            params_i = dynamics.get_mode_dyn_params(u, constants, m_i);
            [mu_x_i, ~]       = filters.ekf_time_update(mu_H_pos(:,m_i+1), cov_H_pos(:,:,m_i+1), params_i);
            params_j = dynamics.get_mode_dyn_params(u, constants, m_j);
            [mu_x_j, cov_x_j] = filters.ekf_time_update(mu_H_pos(:,m_j+1), cov_H_pos(:,:,m_j+1), params_j);
            mu_y_i = params_i.H*mu_x_i;
            mu_y_j = params_j.H*mu_x_j;
            cov_y_j = params_j.H*cov_x_j*params_j.H' + params_j.R;
            distance = get_Mah_distance(mu_y_i, mu_y_j, cov_y_j);
            h_i = Hypothesis(m_i+1); h_j = Hypothesis(m_j+1);
            % distance = distance*sqrt(h_i)*sqrt(h_j);
            metric = metric + distance;
            min_dist = min(min_dist, distance);
            geo_mean = geo_mean*distance;
        end
    end
end
% metric = metric * min_dist; % to guarantee exciting control
metric = geo_mean^(1/(n_Hypo^2));
end

function d = get_Mah_distance(x, mu, cov)
d = sqrt((x-mu)'*cov^(-1)*(x-mu));
end
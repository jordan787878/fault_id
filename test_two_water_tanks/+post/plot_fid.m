function plot_fid(data)
    
plot_fid_filter(data);

plot_fid_control_and_belief(data)

end


function plot_fid_filter(data)

log = "(True): " + num2str(data.true_mode) + " (ID): " + num2str(data.id_mode);
disp(log);

nk = size(data.x_hist, 2);
ks = 0:1:(nk-1);
ny     = size(data.y_hist, 1);

% filter results of the true hypothesis: true state, est. state, and covariance bounds.
figure(1)
x     = data.x_hist;
y     = data.y_hist;
nx = size(x, 1);
if(data.true_mode >= 0)
    plot_mode = data.true_mode+1;
else
    plot_mode = 1;
end
x_est = reshape(data.mu_H_pos_hist(:,plot_mode), nx, []);
x_cov = data.cov_H_pos_hist(:,:,plot_mode);
colorList = lines(nx);
for i = 1:nx
    subplot(2, 1, i);
    hold on;

    % Plot true state
    p1 = plot(ks, x(i,:), 'k--', 'LineWidth', 1.5);
    p1.DisplayName = "true state";

    % Plot noisy measurement
    p2 = scatter(ks, y(i,:), 30, "k", Marker="x", LineWidth=1.5);
    p2.DisplayName = "noisy measurement";

    p3 = plot(ks, x_est(i,:),'LineStyle', '-', 'Color', colorList(i,:), 'LineWidth', 2.0);
    p3.DisplayName = "state estimate";

    sigma = sqrt(x_cov(i:nx:end,i)');
    upper_bound = x_est(i,:) + 3 * sigma;
    lower_bound = x_est(i,:) - 3 * sigma;
    % Create x and y outputs for fill
    x_fill = [ks, fliplr(ks)];                   
    y_fill = [upper_bound, fliplr(lower_bound)]; 
    p4 = fill(x_fill, y_fill, colorList(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    p4.DisplayName = "3 sigma bound";
    ylim([min(lower_bound), max(upper_bound)])
end
% Add x and y labels to the first subplot
subplot(2,1,1)
xlabel('time step, k');
ylabel('Channel 1');
legend('show',Location="best");
% Add x and y labels to the second subplot
subplot(2,1,2)
xlabel('time step, k');
ylabel('Channel 2');



figure(2); 
% n_Hypo = size(data.H_hist, 2) - 1; % ignore last "fake" hypo
n_Hypo = size(data.H_hist, 2);
colorList = lines(n_Hypo);
for m = 0:n_Hypo-1
    inno = reshape(data.inno_H_pri_hist(:,m+1), ny, nk);
    inno_whiten = 0.0*inno;
    S    = data.S_H_pri_hist(:,:,m+1);
    for k = 0:(nk-1)
        inno_k = inno(:,k+1);
        S_k    = S(ny*k+1:ny*k+ny,:);
        inno_whiten_k = inv(chol(S_k, 'lower'))*inno_k;
        inno_whiten(:,k+1) = inno_whiten_k;
    end

    for j = 1:2
    subplot(2,1,j)
    if(m == data.id_mode && m == data.true_mode)
       h = plot(ks, inno_whiten(j,:), 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (ID, True)";
    elseif(m == data.id_mode)
       sc = scatter(ks, inno_whiten(j,:), 150, 'Marker', 'o', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m) +" (ID)";
    elseif(m == data.true_mode)
       h = plot(ks, inno_whiten(j,:), 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (True)";
    else
       sc = scatter(ks, inno_whiten(j,:), 50, 'Marker', 'x', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m);
    end
    end
end

% Add x and y labels to the first subplot
subplot(2,1,1)
xlabel('time step, k');
ylabel('Whiten dy (Channel 1)');
legend('show');

% Add x and y labels to the second subplot
subplot(2,1,2)
xlabel('time step, k');
ylabel('Whiten dY (Channel 2)');

sgtitle(log)

end


function plot_fid_control_and_belief(data)
n_Hypo = size(data.H_hist, 2) - 1; % ignore last "fake" hypo
n_Hypo = size(data.H_hist, 2);
colorList = lines(n_Hypo);
nk = size(data.x_hist, 2);
ks = 0:1:(nk-1);
ny     = size(data.y_hist, 1);

figure(3); 
for m = 0:n_Hypo-1
    belief = data.H_hist(:, m+1);
    subplot(2,1,1)
    if(m == data.id_mode && m == data.true_mode)
       h = plot(ks, belief, 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (ID, True)";
    elseif(m == data.id_mode)
       sc = scatter(ks, belief, 150, 'Marker', 'o', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m) +" (ID)";
    elseif(m == data.true_mode)
       h = plot(ks, belief, 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (True)";
    else
       sc = scatter(ks, belief, 50, 'Marker', 'x', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m);
    end

    subplot(2,1,2)
    plot(ks, data.u_hist, "k-", 'LineWidth', 1.5);
    ylim([0.0, 0.005]); % [tmp hard code]
end

% Add x and y labels to the first subplot
subplot(2,1,1)
xlabel('time step, k');
ylabel('Hypothesis Belief');
legend('show',Location="northwest");
 
% Add x and y labels to the second subplot
subplot(2,1,2)
xlabel('time step, k');
ylabel('control');

log = "(True): " + num2str(data.true_mode) + " (ID): " + num2str(data.id_mode);
sgtitle(log)
end
function plot_fid(outputs)
post.setPublicationDefaults
% % animation
% post.animate(fid);

% logging
disp("True mode: " + num2str(outputs.true_mode) + ", ID mode: " + num2str(outputs.id_mode));

% filter results of the true hypothesis: true state, est. state, and covariance bounds.
figure(2)
t = outputs.t_hist;
t = t(1:outputs.k_end);
x     = outputs.x_hist;
x = x(:,1:outputs.k_end);
y     = outputs.y_hist;
y = y(:,1:outputs.k_end);
nx = size(x, 1);
if(outputs.true_mode >= 0)
    plot_mode = outputs.true_mode+1;
else
    plot_mode = 1;
end
x_est = reshape(outputs.mu_H_pos_hist(:,plot_mode), nx, []);
x_cov = outputs.cov_H_pos_hist(:,:,plot_mode);
colorList = lines(nx);
subplot_order = [1, 3, 5, 2, 4, 6];
for i = 1:nx
    subplot(3, 2, subplot_order(i));
    hold on;
    
    % Plot true state
    p1 = plot(t, x(i,:), 'k--', 'LineWidth', 1.5);
    p1.DisplayName = "true state";

    % Plot noisy measurement
    p2 = scatter(t, y(i,:), 30, "k", Marker="x", LineWidth=1.5);
    p2.DisplayName = "noisy measurement";
    
    % Plot filter state estimate (X_ref + x_est)
    p3 = plot(t, x_est(i,:),'LineStyle', '-', 'Color', colorList(i,:), 'LineWidth', 2.0);
    p3.DisplayName = "state estimate";

    sigma = sqrt(x_cov(i:nx:end,i)');
    upper_bound = x_est(i,:) + 6 * sigma;
    lower_bound = x_est(i,:) - 6 * sigma;
    
    % Create x and y outputs for fill
    x_fill = [t, fliplr(t)];                   
    y_fill = [upper_bound, fliplr(lower_bound)]; 
    
    % Fill the area between the bounds with a translucent patch
    p4 = fill(x_fill, y_fill, colorList(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    p4.DisplayName = "6 sigma bound";

    ylim([min(lower_bound), max(upper_bound)])
    xlabel('t, sec')
    ylabel("x "+num2str(i))
end
% Add x and y labels to the first subplot
subplot(3,2,1)
% xlabel('time, sec');
% ylabel('state (Channel 1)');
legend('show',Location="best");


% filter whiten y inno.
figure(3);
n_Hypo = size(outputs.mu_H_pos_hist,2);
colorList = lines(n_Hypo);
nk = outputs.k_end;
ny     = size(outputs.y_hist, 1);
markerstyle = ["x","s","d","o"];
for m = 0:n_Hypo-1
    inno = reshape(outputs.inno_H_pri_hist(:,m+1), ny, nk);
    inno_whiten = 0.0*inno;
    S    = outputs.S_H_pri_hist(:,:,m+1);
    for k = 0:(nk-1)
        inno_k = inno(:,k+1);
        S_k    = S(ny*k+1:ny*k+ny,:);
        inno_whiten_k = inv(chol(S_k, 'lower'))*inno_k;
        inno_whiten(:,k+1) = inno_whiten_k;
    end
    for j = 1:nx
        subplot(3, 2, subplot_order(j));
        sc = scatter(t, inno_whiten(j,:), 30, 'Marker', 'x', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.0); hold on;
        sc.DisplayName = num2str(m);
        if(m == outputs.true_mode)
            sc.Marker = 'o';
            sc.SizeData = 100;
        end
        xlabel('t, sec')
        ylabel("dy "+num2str(j))
    end
end
% Add x and y labels to the first subplot
subplot(3,2,1)
% xlabel('time, sec');
% ylabel('whiten inno (Channel 1)');
legend('show',Location="best");

% filter hypothesis belief and applied control
figure(4); 
for m = 0:n_Hypo-1
    belief = outputs.H_hist(:, m+1);
    subplot(2,1,1)
    if(m == outputs.id_mode && m == outputs.true_mode)
       h = plot(t, belief, 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (ID, True)";
    elseif(m == outputs.id_mode)
       sc = scatter(t, belief, 500, 'Marker', 'o', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m) +" (ID)";
    elseif(m == outputs.true_mode)
       h = plot(t, belief, 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (True)";
    else
       sc = scatter(t, belief, 150, 'Marker', 'x', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m);
    end

    subplot(2,1,2)
    u_hist = outputs.u_hist(:,1:outputs.k_end);
    plot(t, u_hist, 'LineWidth', 1.5);
end

% Add x and y labels to the first subplot
subplot(2,1,1)
xlabel('t, sec');
ylabel('Hypothesis Belief');
legend('show',Location="northwest");
 
% Add x and y labels to the second subplot
subplot(2,1,2)
xlabel('t, sec');
ylabel('control');
legend('Channel 1', 'Channel 2', 'Channel 3', Location="Best")

end

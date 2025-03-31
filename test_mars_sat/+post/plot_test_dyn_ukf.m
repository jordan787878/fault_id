function plot_test_dyn_ukf(fid)

% post.animate(fid);

figure(2)
t = fid.outputs.t_hist;
x = fid.outputs.x_hist;
y = fid.outputs.y_hist;
u = fid.outputs.u_hist;
nx = size(x, 1);
m = 0;
x_est = reshape(fid.outputs.mu_H_pos_hist(:,m+1), nx, []);
x_cov = fid.outputs.cov_H_pos_hist(:,:,m+1);
colorList = lines(nx);
subplot_order = [1, 3, 5, 2, 4, 6];
for i = 1:nx
    subplot(3, 2, subplot_order(i));
    hold on;
    
    % Plot true state
    plot(t, x(i,:), 'k', 'LineWidth', 1.5);
    
    % Plot filter state estimate (X_ref + x_est)
    plot(t, x_est(i,:),'LineStyle', '--', 'Color', colorList(i,:), 'LineWidth', 2.0);

    scatter(t(1:10:end), y(i,1:10:end), 20, 'x', 'MarkerEdgeColor', colorList(i,:), 'MarkerFaceColor', colorList(i,:))

    sigma = sqrt(x_cov(i:nx:end,i)');
    upper_bound = x_est(i,:) + 6 * sigma;
    lower_bound = x_est(i,:) - 6 * sigma;
    
    % Create x and y data for fill
    x_fill = [t, fliplr(t)];                   
    y_fill = [upper_bound, fliplr(lower_bound)]; 
    
    % Fill the area between the bounds with a translucent patch
    fill(x_fill, y_fill, colorList(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    xlabel("t, sec")
    ylabel("x"+num2str(i))
end


figure(3);
data = fid.outputs;
n_Hypo = size(data.mu_H_pos_hist,2);
colorList = lines(n_Hypo);
nk = length(t);
ny     = size(data.y_hist, 1);
m = 0;
inno = reshape(data.inno_H_pri_hist(:,m+1), ny, nk);
inno_whiten = 0.0*inno;
S    = data.S_H_pri_hist(:,:,m+1);
for k = 0:(nk-1)
    inno_k = inno(:,k+1);
    S_k    = S(ny*k+1:ny*k+ny,:);
    inno_whiten_k = inv(chol(S_k, 'lower'))*inno_k;
    inno_whiten(:,k+1) = inno_whiten_k;
end
for j = 1:nx
    subplot(3, 2, subplot_order(j));
    sc = scatter(t, inno_whiten(j,:), 10, 'Marker', 'x', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.0); hold on;
    sc.DisplayName = num2str(m);
    xlabel("t, sec")
    ylabel("dy "+num2str(j))
end

% filter hypothesis belief and applied control
figure(4); 
plot(t, u, 'LineWidth', 1.5);
xlabel('t, sec');
ylabel('control');
legend('Channel 1', 'Channel 2', 'Channel 3', Location="Best")


end

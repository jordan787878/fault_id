function plot_dataset1b1()
    post.setPublicationDefaults
    
    datasetFolder = "data/dataset1b1";
    N = 25;
    idx = 3;
    folder_hig = fullfile(datasetFolder, 'rlevel_hig');
    activeFileName  = sprintf('data_act_win%d.mat', N);
    passiveFileName = sprintf('data_pas_win%d.mat', N);
    
    pas = load(fullfile(folder_hig, passiveFileName));
    fn = fieldnames(pas);
    pas_outputs = pas.(fn{1}).sim.outputs{idx};
    
    act = load(fullfile(folder_hig, activeFileName));
    fn = fieldnames(act);
    act_outputs = act.(fn{1}).sim.outputs{idx};
    
    % figure(1) control compare
    figure(1)
    color_blue  = [0 0.4470 0.7410];    % MATLAB default blue
    color_red = [0.8500 0.3250 0.0980]; % MATLAB default red
    color_green = [0, 0.5, 0];
    colorList = [color_blue; color_red; color_green];
    global_min = min([pas_outputs.u_hist(:); act_outputs.u_hist(:)]);
    global_max = max([pas_outputs.u_hist(:); act_outputs.u_hist(:)]);
    for i = 1:3
        subplot(3,1,i)
        plot(act_outputs.t_hist, act_outputs.u_hist(i,:), 'LineStyle', '--', Color=colorList(i,:)); hold on;
        plot(pas_outputs.t_hist, pas_outputs.u_hist(i,:), 'Color', colorList(i,:)); hold on;
        
        ylabel("u"+num2str(i))
        % xline(act_outputs.k_end + act_outputs.t_hist(1), 'Color', 'k');
        legend('Active','Passive',Location="southeast")
        title("Control Channel "+num2str(i))
        ylim([global_min, global_max])
    end

    xlabel("k, time step")
    
    figure(2)
    t_min = act_outputs.t_hist(1); t_max = t_min + pas_outputs.k_end;
    subplot(2,1,1)
    plot_belief_traj(act_outputs)
    xlim([t_min, t_max])
    legend('show',Location="northwest");
    xlabel("k, time step")
    ylabel("belief")
    title("Active Hypothesis Belief Trajectory")
    subplot(2,1,2)
    plot_belief_traj(pas_outputs)
    xlim([t_min, t_max])
    legend('show',Location="northwest");
    xlabel("k, time step")
    ylabel("belief")
    title("Passive Hypothesis Belief Trajectory")

    figure(3)
    gap = 100;
    h_act = visualizeAttitude(act_outputs.bn_hist, act_outputs.rlmo_hist, gap, 3, 0.5, ...
                               {'Active Body X','Active Body Y','Active Body Z'}, '--');
    hold on
    h_pas = visualizeAttitude(pas_outputs.bn_hist, pas_outputs.rlmo_hist, gap, 1, 3.5, ...
                               {'Passive Body X','Passive Body Y','Passive Body Z'}, ':');
    
    % Now combine handles. For example, if you want one legend per dataset:
    legend([h_act(1), h_pas(1)], {'Active Attitude Body X', 'Passive Attitude Body X'}, 'Location', 'best');
    view(25, 11)

    % Ensure that the figure is not cropped when saving.
    set(gcf, 'PaperPositionMode', 'auto');
    
    % Save the figure as a PDF with 300 dpi.
    print(gcf, 'figs/3datt_traj.pdf', '-dpdf', '-r300');

end


%% 
function plot_belief_traj(outputs)
t = outputs.t_hist;
t = t(1:outputs.k_end);
n_Hypo = size(outputs.H_hist,2);
colorList = lines(n_Hypo);
for m = 0:n_Hypo-1
    belief = outputs.H_hist(:, m+1);
    if(m == outputs.id_mode && m == outputs.true_mode)
       h = plot(t, belief, 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (ID, True)";
    elseif(m == outputs.id_mode)
       sc = scatter(t, belief, 150, 'Marker', 'o', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m) +" (ID)";
    elseif(m == outputs.true_mode)
       h = plot(t, belief, 'LineStyle', '--', 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', colorList(m+1, :)); hold on;
       h.DisplayName = num2str(m) +" (True)";
    else
       sc = scatter(t, belief, 50, 'Marker', 'x', 'MarkerEdgeColor', colorList(m+1,:), 'LineWidth', 1.5); hold on;
       sc.DisplayName = num2str(m);
    end
end
end

function h = visualizeAttitude(bn_list, rlmo_list, gap, factor, lw, legendLabels, quiverLineStyle)
% visualizeAttitude Visualizes attitude changes and positions over time.
%
%   visualizeAttitude(bn_list, rlmo_list, gap, factor, lw, legendLabels, quiverLineStyle)
%   creates a 3D plot that shows:
%       - The trajectory of the center positions given by rlmo_list.
%       - Attitude frames (DCMs) from bn_list plotted as body-axis arrows.
%
%   Inputs:
%       bn_list       - 3x3xN array where bn_list(:,:,i) is the DCM at time i.
%       rlmo_list     - 3xN array where each column is the (x,y,z) position.
%       gap           - Sampling gap; e.g., gap = 2 plots every other attitude frame
%                       (i.e., indices 1:gap:(N-1), so that for an array of length 21 and
%                       gap=10, you get indices [1, 10, 20]).
%       factor        - A scaling factor applied to the arrow length.
%       lw            - Line width used for the quiver arrows.
%       legendLabels  - A cell array of three strings for the legend, e.g.,
%                       {'Body X-axis','Body Y-axis','Body Z-axis'}.
%       quiverLineStyle - A string specifying the linestyle for the quiver arrows,
%                         e.g. '--', ':', '-.'.
%
%   Example:
%       visualizeAttitude(bn_list, rlmo_list, 2, 1, 2, ...
%                         {'X-axis','Y-axis','Z-axis'}, '--');
%
%   Author: Your Name
%   Date:   Today's Date

% Validate dimensions.
if size(bn_list,1) ~= 3 || size(bn_list,2) ~= 3
    error('bn_list must be a 3x3xN array.');
end
if size(rlmo_list,1) ~= 3
    error('rlmo_list must be a 3xN array.');
end

N = size(bn_list,3);
if size(rlmo_list,2) ~= N
    error('The number of attitude matrices in bn_list must match the number of positions in rlmo_list.');
end

% Compute uniform indices: use indices = [1, gap, 2*gap, ...] provided they are strictly less than N.
if gap > 1
    indices = 1:gap:(N-1);
else
    indices = 1:N;
end

% Set up the axes (assumes the figure is already created or will be created before calling).
hold on;
grid on;
axis equal;
xlabel('X, km');
ylabel('Y, km');
zlabel('Z, km');
title('Attitude Changes Over Time');

% Define a scaling factor for the body axes so they appear clearly.
max_distance = max(vecnorm(rlmo_list, 2, 1));
axis_scale = 0.1 * max_distance;
if axis_scale == 0
    axis_scale = 1;
end
axis_scale = axis_scale * factor;

scatter3(rlmo_list(1, 1),rlmo_list(2, 1),rlmo_list(3, 1), 150, "black","filled", "o");

% Loop over the selected time indices and draw the attitude frames.
for i = indices
    center = rlmo_list(:, i);
    R = bn_list(:, :, i);
    
    % Draw the body x-axis (red)
    qx = quiver3(center(1), center(2), center(3), ...
            axis_scale * R(1,1), axis_scale * R(2,1), axis_scale * R(3,1), ...
            'r', 'LineWidth', lw, 'MaxHeadSize', 0.5);
    % Set linestyle for x-axis.
    set(findobj(qx, 'Type', 'line'), 'LineStyle', quiverLineStyle);
    
    % Draw the body y-axis (green)
    qy = quiver3(center(1), center(2), center(3), ...
            axis_scale * R(1,2), axis_scale * R(2,2), axis_scale * R(3,2), ...
            'g', 'LineWidth', lw, 'MaxHeadSize', 0.5);
    set(findobj(qy, 'Type', 'line'), 'LineStyle', quiverLineStyle);
    
    % Draw the body z-axis (blue)
    qz = quiver3(center(1), center(2), center(3), ...
            axis_scale * R(1,3), axis_scale * R(2,3), axis_scale * R(3,3), ...
            'b', 'LineWidth', lw, 'MaxHeadSize', 0.5);
    set(findobj(qz, 'Type', 'line'), 'LineStyle', quiverLineStyle);
end

% Create dummy handles for the legend.
h1 = plot3(NaN, NaN, NaN, 'r', 'LineWidth', lw);
h2 = plot3(NaN, NaN, NaN, 'g', 'LineWidth', lw);
h3 = plot3(NaN, NaN, NaN, 'b', 'LineWidth', lw);
h = [h1, h2, h3];
axis equal;
end

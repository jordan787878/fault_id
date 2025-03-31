function animate_twosat(act_outputs, pas_outputs)
% animate_positions_two_satellites Animates attitude evolution for two satellites.
%
%   animate_positions_two_satellites(act_outputs, pas_outputs) creates a 3D animation
%   of two satellite outputs (active and passive) on the same figure. The attitude sets
%   are distinguished by different line styles and line widths.
%
%   Each input is a structure with fields:
%       rlmo_hist   - 3xN array of the center position history.
%       rgmo_hist   - 3xN array of secondary positions.
%       bn_hist     - 3x3xN array of attitude DCMs.
%       t_hist      - 1xN array of time stamps.
%
%   Author: Your Name
%   Date: Today's Date

%% Extract data for active satellite
R1_act = act_outputs.rlmo_hist;
R2_act = act_outputs.rgmo_hist;
BN_act = act_outputs.bn_hist;
T_act  = act_outputs.t_hist;

%% Extract data for passive satellite
R1_pas = pas_outputs.rlmo_hist;
R2_pas = pas_outputs.rgmo_hist;
BN_pas = pas_outputs.bn_hist;
T_pas  = pas_outputs.t_hist;

%% Determine number of frames (use minimum length if they differ)
N_act = size(R1_act,2);
N_pas = size(R1_pas,2);
N = min(N_act, N_pas);

%% Create VideoWriter object
writerObj = VideoWriter("figs/anim_two_sat.mp4", 'MPEG-4');
writerObj.FrameRate = 24;
open(writerObj);

%% Create figure and set up axes
fig = figure('Position', [100, 100, 800, 600]);
hold on;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Animation of Active and Passive Satellites');

% Compute axis limits from union of positions from both satellites
margin = 3000;
allX = [R1_act(1,:) R2_act(1,:) R1_pas(1,:) R2_pas(1,:)];
allY = [R1_act(2,:) R2_act(2,:) R1_pas(2,:) R2_pas(2,:)];
allZ = [R1_act(3,:) R2_act(3,:) R1_pas(3,:) R2_pas(3,:)];
x_min = min(allX)-margin;  x_max = max(allX)+margin;
y_min = min(allY)-margin;  y_max = max(allY)+margin;
z_min = min(allZ)-margin;  z_max = max(allZ)+margin;
xlim([x_min, x_max]);
ylim([y_min, y_max]);
zlim([z_min, z_max]);

view(105,25);
axis equal;

% txt
textbox_handle = annotation('textbox', [0.1, 0.8, 0.1, 0.1], 'String', '', 'FitBoxToText', 'on');

% Plot the Earth as a filled black circle
scatter3(0, 0, 0, 200, "k", "filled");

%% Plot initial positions
h_R1_act = plot3(R1_act(1,1), R1_act(2,1), R1_act(3,1), 'bo', 'MarkerSize', 10);
% h_R1_pas = plot3(R1_pas(1,1), R1_pas(2,1), R1_pas(3,1), 'ro', 'MarkerSize', 10);
h_R2 = plot3(R2_pas(1,1), R2_pas(2,1), R2_pas(3,1), 'co', 'MarkerSize', 15, 'LineWidth', 4);

%% Plot initial attitude sets
% Define arrow length for attitude visualization
arrow_length = 3000;
fac = 0.5;
% Active satellite: use solid lines with thicker width
BN_act_1 = plot3([R1_act(1,1), R1_act(1,1)+arrow_length*BN_act(1,1,1)], ...
                 [R1_act(2,1), R1_act(2,1)+arrow_length*BN_act(1,2,1)], ...
                 [R1_act(3,1), R1_act(3,1)+arrow_length*BN_act(1,3,1)], 'r--', 'LineWidth', 1.5);
BN_act_2 = plot3([R1_act(1,1), R1_act(1,1)+arrow_length*BN_act(2,1,1)], ...
                 [R1_act(2,1), R1_act(2,1)+arrow_length*BN_act(2,2,1)], ...
                 [R1_act(3,1), R1_act(3,1)+arrow_length*BN_act(2,3,1)], 'g--', 'LineWidth', 1.5);
BN_act_3 = plot3([R1_act(1,1), R1_act(1,1)+arrow_length*BN_act(3,1,1)], ...
                 [R1_act(2,1), R1_act(2,1)+arrow_length*BN_act(3,2,1)], ...
                 [R1_act(3,1), R1_act(3,1)+arrow_length*BN_act(3,3,1)], 'b--', 'LineWidth', 1.5);

% Passive satellite: use dashed lines with a thicker width for distinction
BN_pas_1 = plot3([R1_pas(1,1), R1_pas(1,1)+fac*arrow_length*BN_pas(1,1,1)], ...
                 [R1_pas(2,1), R1_pas(2,1)+fac*arrow_length*BN_pas(1,2,1)], ...
                 [R1_pas(3,1), R1_pas(3,1)+fac*arrow_length*BN_pas(1,3,1)], 'r', 'LineWidth', 3);
BN_pas_2 = plot3([R1_pas(1,1), R1_pas(1,1)+fac*arrow_length*BN_pas(2,1,1)], ...
                 [R1_pas(2,1), R1_pas(2,1)+fac*arrow_length*BN_pas(2,2,1)], ...
                 [R1_pas(3,1), R1_pas(3,1)+fac*arrow_length*BN_pas(2,3,1)], 'g', 'LineWidth', 3);
BN_pas_3 = plot3([R1_pas(1,1), R1_pas(1,1)+fac*arrow_length*BN_pas(3,1,1)], ...
                 [R1_pas(2,1), R1_pas(2,1)+fac*arrow_length*BN_pas(3,2,1)], ...
                 [R1_pas(3,1), R1_pas(3,1)+fac*arrow_length*BN_pas(3,3,1)], 'b', 'LineWidth', 3);

%% Plot nadir and GMO lines for both satellites
h_nadir_act = plot3([R1_act(1,1), 0], [R1_act(2,1), 0], [R1_act(3,1), 0], 'k:', 'LineWidth', 0.5);
h_gmo_act   = plot3([R1_act(1,1), R2_act(1,1)], [R1_act(2,1), R2_act(2,1)], [R1_act(3,1), R2_act(3,1)], 'c:', 'LineWidth', 1);
h_nadir_pas = plot3([R1_pas(1,1), 0], [R1_pas(2,1), 0], [R1_pas(3,1), 0], 'k:', 'LineWidth', 0.5);
h_gmo_pas   = plot3([R1_pas(1,1), R2_pas(1,1)], [R1_pas(2,1), R2_pas(2,1)], [R1_pas(3,1), R2_pas(3,1)], 'c:', 'LineWidth', 1);

%% Animation Loop
frame_gap = 10;
for i = 2:frame_gap:N
    % Update active satellite positions and attitude
    set(h_R1_act, 'XData', R1_act(1,i), 'YData', R1_act(2,i), 'ZData', R1_act(3,i));
    set(BN_act_1, 'XData', [R1_act(1,i), R1_act(1,i)+arrow_length*BN_act(1,1,i)], ...
                  'YData', [R1_act(2,i), R1_act(2,i)+arrow_length*BN_act(1,2,i)], ...
                  'ZData', [R1_act(3,i), R1_act(3,i)+arrow_length*BN_act(1,3,i)]);
    set(BN_act_2, 'XData', [R1_act(1,i), R1_act(1,i)+arrow_length*BN_act(2,1,i)], ...
                  'YData', [R1_act(2,i), R1_act(2,i)+arrow_length*BN_act(2,2,i)], ...
                  'ZData', [R1_act(3,i), R1_act(3,i)+arrow_length*BN_act(2,3,i)]);
    set(BN_act_3, 'XData', [R1_act(1,i), R1_act(1,i)+arrow_length*BN_act(3,1,i)], ...
                  'YData', [R1_act(2,i), R1_act(2,i)+arrow_length*BN_act(3,2,i)], ...
                  'ZData', [R1_act(3,i), R1_act(3,i)+arrow_length*BN_act(3,3,i)]);
    set(h_nadir_act, 'XData', [R1_act(1,i), 0], 'YData', [R1_act(2,i), 0], 'ZData', [R1_act(3,i), 0]);
    set(h_gmo_act, 'XData', [R1_act(1,i), R2_act(1,i)], 'YData', [R1_act(2,i), R2_act(2,i)], 'ZData', [R1_act(3,i), R2_act(3,i)]);
    
    % Update passive satellite positions and attitude
    % set(h_R1_pas, 'XData', R1_pas(1,i), 'YData', R1_pas(2,i), 'ZData', R1_pas(3,i));
    set(h_R2, 'XData', R2_pas(1,i), 'YData', R2_pas(2,i), 'ZData', R2_pas(3,i));
    set(BN_pas_1, 'XData', [R1_pas(1,i), R1_pas(1,i)+fac*arrow_length*BN_pas(1,1,i)], ...
                  'YData', [R1_pas(2,i), R1_pas(2,i)+fac*arrow_length*BN_pas(1,2,i)], ...
                  'ZData', [R1_pas(3,i), R1_pas(3,i)+fac*arrow_length*BN_pas(1,3,i)]);
    set(BN_pas_2, 'XData', [R1_pas(1,i), R1_pas(1,i)+fac*arrow_length*BN_pas(2,1,i)], ...
                  'YData', [R1_pas(2,i), R1_pas(2,i)+fac*arrow_length*BN_pas(2,2,i)], ...
                  'ZData', [R1_pas(3,i), R1_pas(3,i)+fac*arrow_length*BN_pas(2,3,i)]);
    set(BN_pas_3, 'XData', [R1_pas(1,i), R1_pas(1,i)+fac*arrow_length*BN_pas(3,1,i)], ...
                  'YData', [R1_pas(2,i), R1_pas(2,i)+fac*arrow_length*BN_pas(3,2,i)], ...
                  'ZData', [R1_pas(3,i), R1_pas(3,i)+fac*arrow_length*BN_pas(3,3,i)]);
    set(h_nadir_pas, 'XData', [R1_pas(1,i), 0], 'YData', [R1_pas(2,i), 0], 'ZData', [R1_pas(3,i), 0]);
    set(h_gmo_pas, 'XData', [R1_pas(1,i), R2_pas(1,i)], 'YData', [R1_pas(2,i), R2_pas(2,i)], 'ZData', [R1_pas(3,i), R2_pas(3,i)]);
    
    % Keep axis limits fixed
    xlim([x_min, x_max]); ylim([y_min, y_max]); zlim([z_min, z_max]);
    time_text = sprintf("T: "+T_pas(i));
    textbox_handle.String = time_text;
    
    % Write current frame to the video
    writeVideo(writerObj, getframe(fig));
    pause(0.001);
    drawnow;
end

close(writerObj);
end

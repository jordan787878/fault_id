function animate_positions(outputs)
R1 = outputs.rlmo_hist;
R2 = outputs.rgmo_hist;
BN = outputs.bn_hist;
mode_data = outputs.ctrlmode_hist;
T = outputs.t_hist;

    % Create a VideoWriter object
    writerObj = VideoWriter("figs/anim.mp4", 'MPEG-4');
    writerObj.FrameRate = 24;  % Set the frame rate (adjust as needed)
    open(writerObj);

    % R1 and R2 are 3x6500 matrices representing the position history of objects
    % Create a figure
    fig=figure('Position', [100, 100, 800, 600]);
    hold on;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Animation');
    
    margin = 1000;
    R1_min = min([R1(1,:), R2(1,:)])-margin;
    R1_max = max([R1(1,:), R2(1,:)])+margin;
    R2_min = min([R1(2,:), R2(2,:)])-margin;
    R2_max = max([R1(2,:), R2(2,:)])+margin;
    R3_min = min([R1(3,:), R2(3,:)])-margin;
    R3_max = max([R1(3,:), R2(3,:)])+margin;
    % Set axis limits based on the position data
    xlim([R1_min, R1_max])
    ylim([R2_min, R2_max])
    zlim([R3_min, R3_max])

    view(105, 25)
    axis equal;

    % Plot the earth
    scatter3(0,0,0,200,"k","filled")

    % txt
    textbox_handle = annotation('textbox', [0.1, 0.8, 0.1, 0.1], 'String', '', 'FitBoxToText', 'on');

    % Plot initial positions
    h_R1 = plot3(R1(1,1), R1(2,1), R1(3,1), 'ro', 'MarkerSize', 10);
    h_R2 = plot3(R2(1,1), R2(2,1), R2(3,1), 'bo', 'MarkerSize', 10);
    % legend('R1', 'R2');

    % attitude set (body-x, body-y, body-z)
    length=1000;
    BN_1 = plot3([R1(1,1), R1(1,1)+length*BN(1,1,1)], ...
                 [R1(2,1), R1(2,1)+length*BN(1,2,1)], ...
                 [R1(3,1), R1(3,1)+length*BN(1,3,1)], 'r', 'LineWidth',2);
    BN_2 = plot3([R1(1,1), R1(1,1)+length*BN(2,1,1)], ...
                 [R1(2,1), R1(2,1)+length*BN(2,2,1)], ...
                 [R1(3,1), R1(3,1)+length*BN(2,3,1)], 'g', 'LineWidth',2);
    BN_3 = plot3([R1(1,1), R1(1,1)+length*BN(3,1,1)], ...
                 [R1(2,1), R1(2,1)+length*BN(3,2,1)], ...
                 [R1(3,1), R1(3,1)+length*BN(3,3,1)], 'b', 'LineWidth',2);

    h_nadir = plot3([R1(1,1), 0],[R1(2,1), 0],[R1(3,1), 0], 'r--', 'LineWidth',0.5);
    h_gmo = plot3([R1(1,1), R2(1,1)],[R1(2,1), R2(2,1)],[R1(3,1), R2(3,1)], 'r--', 'LineWidth',0.5);

    % Animation loop
    gap = 10;
    for i = 2:gap:size(R1, 2)
        % Update textbox with mode data
        if (mode_data(i) == 0)
            m = "Sun-Point";
        elseif(mode_data(i) == 2)
            m = "Nadir-Point";
        else
            m = "GMO-Point";
        end

        % Update positions
        set(h_R1, 'XData', R1(1,i), 'YData', R1(2,i), 'ZData', R1(3,i));
        set(h_R2, 'XData', R2(1,i), 'YData', R2(2,i), 'ZData', R2(3,i));
        set(BN_1, 'XData', [R1(1,i), R1(1,i)+length*BN(1,1,i)], ...
                  'YData', [R1(2,i), R1(2,i)+length*BN(1,2,i)], ...
                  'ZData', [R1(3,i), R1(3,i)+length*BN(1,3,i)]);
        set(BN_2, 'XData', [R1(1,i), R1(1,i)+length*BN(2,1,i)], ...
                  'YData', [R1(2,i), R1(2,i)+length*BN(2,2,i)], ...
                  'ZData', [R1(3,i), R1(3,i)+length*BN(2,3,i)]);
        set(BN_3, 'XData', [R1(1,i), R1(1,i)+length*BN(3,1,i)], ...
                  'YData', [R1(2,i), R1(2,i)+length*BN(3,2,i)], ...
                  'ZData', [R1(3,i), R1(3,i)+length*BN(3,3,i)]);
        set(h_nadir, 'XData', [R1(1,i), 0], ...
                     'YData', [R1(2,i), 0], ...
                     'ZData', [R1(3,i), 0]);
        set(h_gmo, 'XData', [R1(1,i), R2(1,i)], ...
                     'YData', [R1(2,i), R2(2,i)], ...
                     'ZData', [R1(3,i), R2(3,i)]);
        
        axis equal;
        xlim([R1_min, R1_max])
        ylim([R2_min, R2_max])
        zlim([R3_min, R3_max])

        % Update textbox with mode data
        mode_text = sprintf("T: "+T(i)+", Mode: "+m);
        textbox_handle.String = mode_text;
        % if (i == 2+10*50 || i == 2+10*210 || i == 2+10*350 || i==2+10*450)
        %     filename = sprintf('figs/frame_%04d.png', i);
        %     saveas(gcf, filename);
        % end

        % Write the current frame to the video
        writeVideo(writerObj, getframe(fig));

        % Pause to control animation speed (adjust as needed)
        pause(0.001);
        drawnow;
    end

    % Close the VideoWriter object
    close(writerObj);
end
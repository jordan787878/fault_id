function plot_dataset(datasetFolder)
post.setPublicationDefaults
    % plot_dataset calls the modular plotting function for two different valueFields.
    %
    % INPUT:
    %   datasetFolder - (string) the dataset folder path, e.g., 'data/dataset1',
    %                   'data/dataset2', or 'data/dataset3'
    
    % Plot for 'fail_percent'
    figure(1);
    plot_dataset_modular(datasetFolder, 'fail_percent', 'Fail Percentage %', 'Fail Percentage vs Moving Window & (Noise Level)');
    
    % Plot for 'success_id_time_for_known_fault'
    figure(2);
    plot_dataset_modular(datasetFolder, 'success_id_time_for_known_fault', 'Average ID Time Duration', 'Average ID Time vs Moving Window & (Noise Level)');
end

function plot_dataset_modular(datasetFolder, valueField, ylabelStr, titleStr)
% plot_dataset_modular plots a chosen value versus Moving Window N for three noise levels.
%
% USAGE:
%   plot_dataset_modular(datasetFolder, valueField, ylabelStr, titleStr)
%
% INPUTS:
%   datasetFolder - (string) the base folder of the dataset, e.g., 'data/dataset1'
%                   (the function expects subfolders 'rlevel_low', 'rlevel_mid', and 'rlevel_hig' inside this folder)
%   valueField    - (string) the field name inside sim.(metaField) to plot, e.g.,
%                   'fail_percent' or 'success_id_time_for_known_fault'
%   ylabelStr     - (string) label for the y-axis.
%   titleStr      - (string) the title for the plot.
%
% The function loads data from the three folders:
%   [datasetFolder]/rlevel_low
%   [datasetFolder]/rlevel_mid
%   [datasetFolder]/rlevel_hig
%
% For each noise level, it extracts the specified value from both active and passive datasets
% for Moving Window N values. For datasetFolder 'data/dataset3', additional N values of 15, 20, and 25
% are used.

    metaField = 'outputs_meta';

    % Define folder paths for the three noise levels using the datasetFolder input
    folder_low = fullfile(datasetFolder, 'rlevel_low');
    folder_mid = fullfile(datasetFolder, 'rlevel_mid');
    folder_hig = fullfile(datasetFolder, 'rlevel_hig');

    % Determine the vector of window values based on datasetFolder
    N_vals = [1:12, 15:5:40];
    N_data = length(N_vals);

    % Preallocate arrays for the chosen values
    active_val_low  = zeros(1, N_data);
    passive_val_low = zeros(1, N_data);
    active_val_mid  = zeros(1, N_data);
    passive_val_mid = zeros(1, N_data);
    active_val_hig  = zeros(1, N_data);
    passive_val_hig = zeros(1, N_data);

    act_violate = [];
    pas_violate = [];

    % Loop over the index (i) corresponding to each N value
    for i = 1:N_data
        N = N_vals(i);
        % Construct file names for active and passive data
        activeFileName  = sprintf('data_act_win%d.mat', N);
        passiveFileName = sprintf('data_pas_win%d.mat', N);

        % ----- Load rlevel_low Data -----
        % Active
        S_active_low = load(fullfile(folder_low, activeFileName));
        fn = fieldnames(S_active_low);
        dataActive_low = S_active_low.(fn{1});
        active_val_low(i) = dataActive_low.sim.(metaField).(valueField);
        act_violate(end+1)  = dataActive_low.sim.(metaField).violate_percent_total;

        % Passive
        S_passive_low = load(fullfile(folder_low, passiveFileName));
        fn = fieldnames(S_passive_low);
        dataPassive_low = S_passive_low.(fn{1});
        passive_val_low(i) = dataPassive_low.sim.(metaField).(valueField);
        pas_violate(end+1)  = dataPassive_low.sim.(metaField).violate_percent_total;

        % ----- Load rlevel_mid Data -----
        % Active
        S_active_mid = load(fullfile(folder_mid, activeFileName));
        fn = fieldnames(S_active_mid);
        dataActive_mid = S_active_mid.(fn{1});
        active_val_mid(i) = dataActive_mid.sim.(metaField).(valueField);
        act_violate(end+1)  = dataActive_mid.sim.(metaField).violate_percent_total;

        % Passive
        S_passive_mid = load(fullfile(folder_mid, passiveFileName));
        fn = fieldnames(S_passive_mid);
        dataPassive_mid = S_passive_mid.(fn{1});
        passive_val_mid(i) = dataPassive_mid.sim.(metaField).(valueField);
        pas_violate(end+1)  = dataPassive_mid.sim.(metaField).violate_percent_total;

        % ----- Load rlevel_hig Data -----
        % Active
        S_active_hig = load(fullfile(folder_hig, activeFileName));
        fn = fieldnames(S_active_hig);
        dataActive_hig = S_active_hig.(fn{1});
        active_val_hig(i) = dataActive_hig.sim.(metaField).(valueField);
        act_violate(end+1)  = dataActive_hig.sim.(metaField).violate_percent_total;

        % Passive
        S_passive_hig = load(fullfile(folder_hig, passiveFileName));
        fn = fieldnames(S_passive_hig);
        dataPassive_hig = S_passive_hig.(fn{1});
        passive_val_hig(i) = dataPassive_hig.sim.(metaField).(valueField);
        pas_violate(end+1)  = dataPassive_hig.sim.(metaField).violate_percent_total;
    end

    % Display overall constraint violations
    disp("active overall constraint violation " + num2str(mean(act_violate)) + " %")
    disp("passive overall constraint violation " + num2str(mean(pas_violate)) + " %")

    % Define Colors (active: blue, passive: red)
    color_active  = [0 0.4470 0.7410];    % MATLAB default blue
    color_passive = [0.8500 0.3250 0.0980]; % MATLAB default red
    % color_dataset4 = [0, 0.5, 0];

    % Create the figure and fill the area for Active and Passive data
    hold on;

    % --- Active Data: Compute overall min and max across low, mid, and hig ---
    active_lower = min([active_val_low; active_val_mid; active_val_hig], [], 1);
    active_upper = max([active_val_low; active_val_mid; active_val_hig], [], 1);
    x_fill = [N_vals, fliplr(N_vals)];
    y_fill_active = [active_lower, fliplr(active_upper)];
    fill(x_fill, y_fill_active, color_active, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % --- Passive Data: Compute overall min and max across low, mid, and hig ---
    passive_lower = min([passive_val_low; passive_val_mid; passive_val_hig], [], 1);
    passive_upper = max([passive_val_low; passive_val_mid; passive_val_hig], [], 1);
    y_fill_passive = [passive_lower, fliplr(passive_upper)];
    fill(x_fill, y_fill_passive, color_passive, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Overlay the Individual Curves for Clarity
    % Active curves for low, mid, and hig datasets
    plot(N_vals, active_val_low, '--o', 'Color', color_active, 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(N_vals, active_val_mid, '--s', 'Color', color_active, 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(N_vals, active_val_hig, '--^', 'Color', color_active, 'LineWidth', 1.0, 'MarkerSize', 5);

    % Passive curves for low, mid, and hig datasets
    plot(N_vals, passive_val_low, '--o', 'Color', color_passive, 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(N_vals, passive_val_mid, '--s', 'Color', color_passive, 'LineWidth', 1.0, 'MarkerSize', 5);
    plot(N_vals, passive_val_hig, '--^', 'Color', color_passive, 'LineWidth', 1.0, 'MarkerSize', 5);

    % Labeling and Title
    xlabel('Moving Window N');
    ylabel(ylabelStr);
    title(titleStr);
    legend({...
        'Active (Noise Low)','Active (Noise Mid)','Active (Noise High)',...
        'Passive (Noise Low)','Passive (Noise Mid)','Passive (Noise High)'}, 'Location', 'northeast');
    grid on;
    hold off;
end

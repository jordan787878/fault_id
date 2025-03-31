function plot_ablation()
    post.setPublicationDefaults
    figure(1);
    plot_dataset_modular('fail_percent', 'Fail Rate %', 'Fail Rate vs Moving Window (Noise Level)');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 11 9]);  % Adjust dimensions as needed.
    set(gcf, 'PaperSize', [11 9]);         % Match this to PaperPosition.
    print(gcf, 'figs/dataset3b_failpercent_ab.pdf', '-dpdf', '-r300');
end

function plot_dataset_modular(valueField, ylabelStr, titleStr)
    metaField = 'outputs_meta';

    % Determine the vector of window values based on datasetFolder
    N_vals = [1:12, 15:5:40];
    N_data = length(N_vals);
    hold on;

    % Define Colors (active: blue, passive: red)
    color_active  = [0 0.4470 0.7410];    % MATLAB default blue
    color_passive = [0.8500 0.3250 0.0980]; % MATLAB default red
    color_ab1 = [0, 0.5, 0];
    color_ab2 = [0.4940, 0.1840, 0.5560];

    % passive
    for i = 1:N_data
        N = N_vals(i);
        passiveFileName = sprintf('data_pas_win%d.mat', N);
        S = load(fullfile("data/dataset3b/rlevel_hig", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_hig(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset3b/rlevel_mid", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_mid(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset3b/rlevel_low", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_low(i) = data.sim.(metaField).(valueField);
    end
    lower = min([val_low; val_mid; val_hig], [], 1);
    upper = max([val_low; val_mid; val_hig], [], 1);
    x_fill = [N_vals, fliplr(N_vals)];
    y_fill = [lower, fliplr(upper)];
    fill(x_fill, y_fill, color_passive, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(N_vals, val_hig, '--o', 'Color', color_passive, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_mid, '--s', 'Color', color_passive, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_low, '--^', 'Color', color_passive, 'LineWidth', 1.5, 'MarkerSize', 5);

    % active
    for i = 1:N_data
        N = N_vals(i);
        passiveFileName = sprintf('data_act_win%d.mat', N);
        S = load(fullfile("data/dataset3b/rlevel_hig", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_hig(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset3b/rlevel_mid", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_mid(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset3b/rlevel_low", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_low(i) = data.sim.(metaField).(valueField);
    end
    lower = min([val_low; val_mid; val_hig], [], 1);
    upper = max([val_low; val_mid; val_hig], [], 1);
    x_fill = [N_vals, fliplr(N_vals)];
    y_fill = [lower, fliplr(upper)];
    fill(x_fill, y_fill, color_active, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(N_vals, val_hig, '--o', 'Color', color_active, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_mid, '--s', 'Color', color_active, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_low, '--^', 'Color', color_active, 'LineWidth', 1.5, 'MarkerSize', 5);

    % ab1: w/o renormalization
    for i = 1:N_data
        N = N_vals(i);
        passiveFileName = sprintf('data_act_win%d.mat', N);
        S = load(fullfile("data/dataset4/rlevel_hig", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_hig(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset4/rlevel_mid", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_mid(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset4/rlevel_low", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_low(i) = data.sim.(metaField).(valueField);
    end
    lower = min([val_low; val_mid; val_hig], [], 1);
    upper = max([val_low; val_mid; val_hig], [], 1);
    x_fill = [N_vals, fliplr(N_vals)];
    y_fill = [lower, fliplr(upper)];
    fill(x_fill, y_fill, color_ab1, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(N_vals, val_hig, '--o', 'Color', color_ab1, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_mid, '--s', 'Color', color_ab1, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_low, '--^', 'Color', color_ab1, 'LineWidth', 1.5, 'MarkerSize', 5);

    % ab1: w/o renormalization
    for i = 1:N_data
        N = N_vals(i);
        passiveFileName = sprintf('data_act_win%d.mat', N);
        S = load(fullfile("data/dataset5/rlevel_hig", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_hig(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset5/rlevel_mid", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_mid(i) = data.sim.(metaField).(valueField);
        S = load(fullfile("data/dataset5/rlevel_low", passiveFileName));
        fn = fieldnames(S);
        data = S.(fn{1});
        val_low(i) = data.sim.(metaField).(valueField);
    end
    lower = min([val_low; val_mid; val_hig], [], 1);
    upper = max([val_low; val_mid; val_hig], [], 1);
    x_fill = [N_vals, fliplr(N_vals)];
    y_fill = [lower, fliplr(upper)];
    fill(x_fill, y_fill, color_ab2, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(N_vals, val_hig, '--o', 'Color', color_ab2, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_mid, '--s', 'Color', color_ab2, 'LineWidth', 1.5, 'MarkerSize', 5);
    plot(N_vals, val_low, '--^', 'Color', color_ab2, 'LineWidth', 1.5, 'MarkerSize', 5);

    % Labeling and Title
    xlabel('Moving Window N');
    ylabel(ylabelStr);
    title(titleStr);
    legend({...
        'Passive (hig)','Passive (mid)','Passive (low)', ...
        'Active (hig)','Active (mid)','Active (low)', ...
        'w/o Renorm (hig)', 'w/o Renorm (mid)', 'w/o Renorm (low)', ...
        'w/o Renorm & Reject (hig)', 'w/o Renorm & Reject (mid)', 'w/o Renorm & Reject (low)'}, ...
        'Location', 'northeast', 'NumColumns', 2, 'FontSize',26);
    grid on;
    hold off;
end

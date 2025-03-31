function main_meta()

% meta script to run main over mulitple trials with different meta_input
% make sure to set the detailed configuraiton in main

meta_inputs_moving_window = [1:12,15:5:40];

% meta_inputs_moving_window = [1,25:25:150]; % only for dataset6

for i =1:length(meta_inputs_moving_window)

    moving_window = meta_inputs_moving_window(i);

    main(moving_window);
end

end
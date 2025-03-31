function main_meta()
% meta script to run main_trials over mulitple trials with different meta_input
% make sure to set the detailed configuraiton in main_trials
meta_inputs_moving_window = [1,5:5:40];
for i =1:length(meta_inputs_moving_window)
    moving_window = meta_inputs_moving_window(i);
    main(moving_window);
end
end
function save_data(sim, fid, targetPath)

data.sim = sim;
data.fid = fid;
if fid.config.action_flag
    status = 'act';
else
    status = 'pas';
end

% Create the variable name string, e.g., 'data_act_win10'
varName = sprintf('data_%s_win%d', status, fid.config.moving_window);
% clear status;

% Now create the variable with that name using eval:
eval([varName ' = data;']);
fileName = sprintf('%s.mat', varName);
fullFileName = fullfile(targetPath, fileName);
save(fullFileName, varName);
% clear data; clear varName; clear fileName; clear fullFileName; clear targetPath;

end
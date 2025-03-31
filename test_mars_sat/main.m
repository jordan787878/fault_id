function main(meta_input)

sim.constants.cov_i = blkdiag(1e-4*eye(3), 1e-6*eye(3));
sim.constants.L_ext = [0;0;0];
sim.constants.I_inertia = diag([10, 5, 7.5]);
sim.constants.L = [zeros(3,3); eye(3)];
sim.constants.Q = 1e-9*eye(3);
sim.constants.H = eye(6);
sim.constants.R = 0.05*blkdiag(1e-4*eye(3), 1e-6*eye(3));
sim.constants.controller_K = (1/6)^2/5;
sim.constants.controller_P = 1/6;

sim.config.N_trials = 2;
sim.config.save_data = false;
sim.config.rlevel = "hig";
sim.config.targetPath = fullfile(pwd, 'data', 'dataset1c', 'rlevel_'+sim.config.rlevel);
sim.config.inc_time = 1;

fid.config.action_flag = false;
fid.config.moving_window = meta_input;
fid.config.N_hypo = 4;
fid.config.crit = 0.95;
fid.config.alpha = 0.05;

% initialization
if(sim.config.rlevel == "hig")
    sim.constants.R = 4.0*sim.constants.R;
elseif(sim.config.rlevel == "mid")
    sim.constants.R = 1.0*sim.constants.R;
else
    sim.constants.R = 0.25*sim.constants.R;
end
sim.outputs = cell(1, sim.config.N_trials);

for i = 1:sim.config.N_trials
    rng(i)

    % random true mode
    true_mode = randi([-1,3]);
    if(true_mode < -1) 
       true_mode = -1; 
    end % ensure true_mode index >= -1

    % random start time
    sim.config.start_time = randi([0,7101]); % based on satellite orbit period ~ 7101 seconds.
    sim.config.end_time = sim.config.start_time + 1000;

    [sim, fid] = run_fid(sim, fid, true_mode);

    % [store results] %
    sim.outputs{i} = fid.outputs;
end

% post-process
sim = post.process_outputs_meta(sim);
disp(sim.outputs_meta)

% plot single trial
% post.plot_fid(sim.outputs{1})
% animation for single trial
% if(sim.config.N_trials == 1)
%     % set saving to false
%     sim.config.save_data = false; post.animate(fid);
% end

% save data
if(sim.config.save_data)
    post.save_data(sim, fid, sim.config.targetPath);
end

end



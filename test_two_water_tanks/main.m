function sim = main(meta_input)

% Test passive vs active fault identification
% action_flag: true or false
% rng(0);

% simulation configuration
sim.config.N_trials = 500;
sim.config.rlevel = "mid";
sim.config.save_data = true;
sim.config.targetPath = fullfile(pwd, 'data', 'dataset4', 'rlevel_'+sim.config.rlevel);
sim.config.verbose = false;
sim.config.use_process_noise = true;   
sim.config.use_measure_noise = true;

% fault identification configuration
fid.config.moving_window = meta_input;
fid.config.action_flag = true;
fid.config.renorm_flag = false;
fid.config.reject_flag = true;
fid.config.filter_options = "ekf";
fid.config.likeli_options = "2"; % "1": raw data; "2": transformed (always use "2" for more numerical stability)
fid.config.N_hypo = 5; %[0,1,2,3,"4"]
fid.config.u_default = 0.0005;
fid.config.crit = 0.95;
fid.config.alpha = 0.05;
fid.config.x1_max = 20.0;
fid.config.x2_min = 1.0;

sim.constants.dt = 5;
sim.constants.A = 1.5e-2;
sim.constants.C1L = 1.2e-4;
sim.constants.C12 = 1.2e-4;
sim.constants.C2L = 1.2e-4;
sim.constants.mu_i = [10.0; 8.0];
sim.constants.cov_i = 0.1*eye(2);
sim.constants.L = eye(2);
sim.constants.H   = eye(2);        % y = Hx;
sim.constants.Q = 0.001*eye(2);
if(sim.config.rlevel == "hig")
    sim.constants.R = 4.0*0.01*eye(2);
elseif(sim.config.rlevel == "mid")
    sim.constants.R = 0.01*eye(2);
else
    sim.constants.R = 0.25*0.01*eye(2);
end

sim.outputs = cell(1, sim.config.N_trials);

for i = 1:sim.config.N_trials
    rng(i);
    true_mode = randi([-1,3]); % true mode = {-1,0,1,2,3}, -1: unknown fault
    if(true_mode < -1) 
       true_mode = -1; 
    end % ensure true_mode index >= -1

    % [running fid] %
    [sim, fid] = run_fid_v2(sim, fid, true_mode);

    % [store results] %
    sim.outputs{i} = fid.outputs;
end

%% post-process
sim = post.process_outputs_meta(sim, fid);
disp(sim.outputs_meta)

%% plots
% post.plot_fid(sim.outputs{1})

%% store all data
if(sim.config.save_data)
    post.save_data(sim, fid, sim.config.targetPath);
end

end


function params = get_mode_dyn_params(u, sim, mode)

constants = sim.constants;
params.L = constants.L;
params.Q = constants.Q;
params.H = constants.H;
params.R = constants.R;

% mode: [0,1,2,3]
% mode: 0 (nominal)
params.RW1 = 1.0;
params.RW2 = 1.0;
params.RW3 = 1.0;

% mode: 1 (fualt 1)
if(mode == 1)
    params.RW1 = 0.2;
end
% mode: 2 (fault 2)
if(mode == 2)
    params.RW2 = 0.2;
end
% mode: 3 (fault 3)
if(mode == 3)
    params.RW3 = 0.2;
end

% (unknown fault)
if(mode == -1)
    params.RW1 = 0.1;
    params.RW2 = 0.5;
    params.RW3 = 0.7;
end

end
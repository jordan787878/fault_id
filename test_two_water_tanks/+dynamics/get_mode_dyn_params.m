function params = get_mode_dyn_params(u, constants, mode)

% mode: [0,1,2,3,4]

params.L = constants.L;
params.Q = constants.Q;
params.H = constants.H;
params.R = constants.R;

% mode: 0 (nominal)
params.u = u;
params.A = constants.A;
params.dt = constants.dt;
params.C1L = 0;
params.C12 = constants.C12;
params.C2L = constants.C2L;

% mode: 1 (fualt 1)
if(mode == 1)
    params.C1L = 1.0*constants.C1L;
end

% mode: 2 (fault 2)
if(mode == 2)
    params.C12 = 0.5*constants.C12;
end
% % mode: 2 (fault 2, test for "similar" modeled dynamics)
% if(mode == 2)
%     params.C1L = 1.05*constants.C1L;
% end

% mode: 3 (fault 3)
if(mode == 3)
    params.C2L = 0.5*5.0*constants.C2L;
end
% mode: 4 (fault 4)
if(mode == 4)
    params.u = 0.5*u;
end

% (unknown fault)
if(mode == -1)
    params.C1L = 0.5*constants.C1L;
    params.C12 = 3.0*constants.C12;
    params.C2L = 1.0*constants.C2L;
end

end
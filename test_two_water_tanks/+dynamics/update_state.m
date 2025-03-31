function x_new = update_state(x, u, w, constants, mode)
params = dynamics.get_mode_dyn_params(u, constants, mode);
u = params.u;
C1L = params.C1L;
C12 = params.C12;
C2L = params.C2L;
dt  = params.dt;
A   = params.A;
q1L = C1L*sqrt(x(1));
q12 = C12*sqrt(x(1)-x(2));
q2L = C2L*sqrt(x(2));
x_new = 0*x;
x_new(1) = x(1) + dt*(u-q1L-q12)/A;
x_new(2) = x(2) + dt*(q12-q2L)/A;
x_new = x_new + constants.L*w;
end
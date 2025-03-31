function x_new = time_update_state(sim, x, u, t1, t2, mode)

% get constants
dt = t2 - t1;
constants = sim.constants;

% get parameters
params = dynamics.get_mode_dyn_params(u, sim, mode);
u(1) = params.RW1 * u(1);
u(2) = params.RW2 * u(2);
u(3) = params.RW3 * u(3);

h1 = dynamics.attitude_dyn(x, t1, u, constants);
h2 = dynamics.attitude_dyn(x+0.5*h1*dt, t1+0.5*dt, u, constants);
h3 = dynamics.attitude_dyn(x+0.5*h2*dt, t1+0.5*dt, u, constants);
h4 = dynamics.attitude_dyn(x+h3*dt, t1+dt, u, constants);
x_new = x + dt*(h1 + 2*h2 + 2*h3 + h4)/6;
% ensure short MRP
[x_new, aux] = satsym.ensure_short_MRP(x_new);

end
function x_new = time_update_state(x, params)
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
end
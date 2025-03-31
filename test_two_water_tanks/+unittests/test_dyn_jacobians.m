function test_dyn_jacobians()
% Define symbolic variables
syms x1 x2 u C1L C12 C2L dt A real

% State vector
x = [x1; x2];

q1L = C1L*sqrt(x1);
q12 = C12*sqrt(x1 - x2);
q2L = C2L*sqrt(x2);

% Define discrete-time dynamics function
f = [x1 + (u - q1L - q12)*dt/A
     x2 + (q12 - q2L)*dt/A];

% Compute the symbolic Jacobian of f with respect to the state x
Jx = simplify(jacobian(f, x));

% Compute the symbolic Jacobian of f with respect to the input u
Ju = simplify(jacobian(f, u));

% Display the discrete-time dynamics and the Jacobians
disp('Discrete-time dynamics:');
disp(f);
disp('Jacobian with respect to state x:');
disp(Jx);
disp('Jacobian with respect to input u:');
disp(Ju);
end
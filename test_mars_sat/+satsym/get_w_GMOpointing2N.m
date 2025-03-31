function w = get_w_GMOpointing2N(t)
dt = 0.1;

dRdt = (satsym.get_GMO_pointing_frame(t+dt) - satsym.get_GMO_pointing_frame(t))/dt;
R = satsym.get_GMO_pointing_frame(t); %[RcN]
W_tilde = dRdt*(-R');
w = [-W_tilde(2,3); W_tilde(1,3); -W_tilde(1,2)];
w = R'*w;
end
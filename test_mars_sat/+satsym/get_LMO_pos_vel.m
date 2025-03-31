function [r,v] = get_LMO_pos_vel(t)
LMO.r = 400 + 3396.10;
LMO.euler_angles_0 = deg2rad([20;30;60]);
LMO.orbit_rate = 0.000884797;
[r, v] = satsym.get_inertial_pos_vel(LMO.r, LMO.euler_angles_0, LMO.orbit_rate, t);
end
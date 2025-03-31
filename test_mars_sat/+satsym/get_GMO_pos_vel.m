function [r,v] = get_GMO_pos_vel(t)
GMO.r = 20424.2;
GMO.euler_angles_0 = deg2rad([0;0;250]);
GMO.orbit_rate = 0.0000709003;
[r, v] = satsym.get_inertial_pos_vel(GMO.r, GMO.euler_angles_0, GMO.orbit_rate, t);
end
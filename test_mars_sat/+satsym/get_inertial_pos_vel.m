function [r_N,v_N] = get_inertial_pos_vel(r, euler_angles_0, orbit_rate, tau)
    e1 = euler_angles_0(1); 
    e2 = euler_angles_0(2); 
    e3 = euler_angles_0(3) + orbit_rate*tau;
    HN = M3(e3)*M1(e2)*M3(e1);
    r_B = [r;0;0];
    r_N = HN'*r_B;
    v_B = r*orbit_rate*[0;1;0];
    v_N = HN'*v_B;
end

function M = M1(t)
M = [1, 0, 0;
     0, cos(t), sin(t);
     0, -sin(t), cos(t)];
end

function M = M3(t)
M = [cos(t), sin(t), 0;
     -sin(t), cos(t), 0;
     0, 0, 1];
end
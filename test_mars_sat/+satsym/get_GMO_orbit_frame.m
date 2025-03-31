function HN = get_GMO_orbit_frame(t)
    euler_angles_0 = deg2rad([0;0;250]);
    orbit_rate = 0.0000709003;
    e1 = euler_angles_0(1); 
    e2 = euler_angles_0(2); 
    e3 = euler_angles_0(3) + orbit_rate*t;
    HN = M3(e3)*M1(e2)*M3(e1);
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
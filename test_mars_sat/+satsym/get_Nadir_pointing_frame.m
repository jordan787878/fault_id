function RnN = get_Nadir_pointing_frame(t)
    RnH = [-1, 0, 0;
            0, 1, 0;
            0, 0,-1];
    RnN = RnH * satsym.get_LMO_orbit_frame(t);
end
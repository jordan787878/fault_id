function w = get_w_Nadir2N(t)
    orbit_rate = 0.000884797;
    HN = satsym.get_LMO_orbit_frame(t);
    w = HN'*[0;0;orbit_rate];
end
function x_init = set_x_init(t_i)
    [r_LMO, ~] = satsym.get_LMO_pos_vel(t_i);
    [r_GMO, ~] = satsym.get_GMO_pos_vel(t_i);
    if(r_LMO(2) > 0)
        mode = 0;
    elseif(satsym.get_angle_deg(r_LMO, r_GMO) < 35)
        mode = 1;
    else
        mode = 2;
    end
    % get control u
    if(mode == 0)
        RN = satsym.get_Sun_pointing_frame(t_i);
        w_R2N_N = satsym.get_w_Sun2N(t_i);
    elseif(mode == 1)
        RN = satsym.get_GMO_pointing_frame(t_i);
        w_R2N_N = satsym.get_w_GMOpointing2N(t_i);
    else
        RN = satsym.get_Nadir_pointing_frame(t_i);
        w_R2N_N = satsym.get_w_Nadir2N(t_i);
    end
    MRP_init = satsym.DCM_to_MRP(RN);
    w_init   = w_R2N_N;
    x_init = [MRP_init; w_init];
end
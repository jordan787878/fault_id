function [ctrl_torque, ctrl_mode, r_LMO, r_GMO, BN, w_BR_Body] = get_control_and_mode(t1, x, mode_force, sim)
%
% ---------------------------------------------------------------------
% 
% Description:
%       compute control torque based on the control mode; 
%       the control mode is based on simple rules to do Sunpointing,
%       Communication with another Sat in GMO, or Nadir pointing
% Outputs:
%       ctrl_torque   - control torque command for attitude pointing
%       ctrl_mode     - control mode {0: SunPoint, 1: Communicate, 2: NadirPoint}
%


[r_LMO, ~] = satsym.get_LMO_pos_vel(t1);
[r_GMO, ~] = satsym.get_GMO_pos_vel(t1);
BN = satsym.get_BN(x(1:3));

if(mode_force == 0)
    if(r_LMO(2) > 0)
        mode = 0;
    elseif(satsym.get_angle_deg(r_LMO, r_GMO) < 35)
        mode = 1;
    else
        mode = 2;
    end
else
    mode = mode_force;
end

% get control u
if(mode == 0)
    RN = satsym.get_Sun_pointing_frame(t1);
    w_R2N_N = satsym.get_w_Sun2N(t1);
elseif(mode == 1)
    RN = satsym.get_GMO_pointing_frame(t1);
    w_R2N_N = satsym.get_w_GMOpointing2N(t1);
else
    RN = satsym.get_Nadir_pointing_frame(t1);
    w_R2N_N = satsym.get_w_Nadir2N(t1);
end
[MRP_BR, w_BR_Body] = satsym.get_attitude_error(t1, x(1:3), x(4:6), RN, w_R2N_N);
u_cmd = -sim.constants.controller_K*MRP_BR ...
        -sim.constants.controller_P*w_BR_Body;
ctrl_torque = u_cmd;
ctrl_mode   = mode;

end
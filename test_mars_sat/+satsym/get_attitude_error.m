function [sigma_B2R, w_B2R_body] = get_attitude_error(t, sigma_B2N, w_B2N_body, RN, w_R2N_N)
% change MRP to DCM
BN = satsym.get_BN(sigma_B2N);
BR = BN*RN';
sigma_B2R = satsym.DCM_to_MRP(BR);
if(norm(sigma_B2R)>1)
    % disp("t switch")
    % disp(t)
    % disp(sigma_B2R)
    sigma_B2R = -sigma_B2R/(norm(sigma_B2R)^2);
    % disp(sigma_B2R)
end
% sigma_N2R = DCM_to_MRP(RN');
% sig_d = sigma_N2R;
% sig_dd = sigma_B2N;
% Num = (1-norm(sig_d)^2)*sig_dd +(1-norm(sig_dd)^2)*sig_d - 2*cross(sig_dd, sig_d);
% Den = 1 + (norm(sig_d)^2)*(norm(sig_dd)^2) - 2*dot(sig_d, sig_dd);
% sigma_B2R = Num/Den;
% if(norm(sigma_B2R)>1)
%     sigma_B2R = -sigma_B2R/(norm(sigma_B2R)^2);
% end
w_B2R_body = w_B2N_body - BN*w_R2N_N;
end
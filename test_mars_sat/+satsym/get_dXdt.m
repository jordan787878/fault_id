function dXdt = get_dXdt(x, t, u, param)
disp("error")
L_ext = param.L_ext;
I_inertia = param.I_inertia;

MRP = x(1:3);
w = x(4:6);
B_MRP = (1-dot(MRP,MRP))*eye(3) + 2*satsym.tilde(MRP) + 2*MRP*MRP';
dMRP_dt = B_MRP*w/4;
dw_dt = inv(I_inertia)*(-satsym.tilde(w)*I_inertia*w + u + L_ext);
dXdt = [dMRP_dt; dw_dt];
end
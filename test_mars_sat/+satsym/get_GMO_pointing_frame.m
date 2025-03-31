function RcN = get_GMO_pointing_frame(t)
LMO.r = 400 + 3396.10;
GMO.r = 20424.2;
n3 = [0;0;1];

r_N_LMO = satsym.get_LMO_orbit_frame(t)'*[LMO.r;0;0];
r_N_GMO = satsym.get_GMO_orbit_frame(t)'*[GMO.r;0;0];
dr = r_N_GMO - r_N_LMO;

Rc1 = -dr/(norm(dr));
Rc2 = cross(dr,n3)/(norm(cross(dr,n3)));
Rc3 = cross(Rc1,Rc2);
RcN = [Rc1'; Rc2'; Rc3'];

end
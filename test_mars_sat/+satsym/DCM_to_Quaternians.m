function beta = DCM_to_Quaternians(C)
b0 = sqrt(0.25*(1+trace(C)));
b1 = sqrt(0.25*(1+2*C(1,1)-trace(C)));
b2 = sqrt(0.25*(1+2*C(2,2)-trace(C)));
b3 = sqrt(0.25*(1+2*C(3,3)-trace(C)));

[~, max_index] = max([b0,b1,b2,b3]);
if(max_index == 1)
    b1 = (C(2,3)-C(3,2))/4/b0;
    b2 = (C(3,1)-C(1,3))/4/b0;
    b3 = (C(1,2)-C(2,1))/4/b0;
end
if(max_index == 2)
    b0 = (C(2,3)-C(3,2))/4/b1;
    b2 = (C(1,2)+C(2,1))/4/b1;
    b3 = (C(3,1)+C(1,3))/4/b1;
end
if(max_index == 3)
    b0 = (C(3,1)-C(1,3))/4/b2;
    b1 = (C(1,2)+C(2,1))/4/b2;
    b3 = (C(2,3)+C(3,2))/4/b2;
end
if(max_index == 4)
    b0 = (C(1,2)-C(2,1))/4/b3;
    b1 = (C(3,1)+C(1,3))/4/b3;
    b2 = (C(2,3)+C(3,2))/4/b3;
end
beta = [b0,b1,b2,b3];
if(b0<0)
    beta = -beta;
end
end
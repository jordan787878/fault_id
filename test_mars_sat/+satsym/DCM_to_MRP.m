function MRP = DCM_to_MRP(DCM)
EP = satsym.DCM_to_Quaternians(DCM);

if(abs(EP(1)+1) < 0.001)
    disp("singular MRP")
end

MRP = [EP(2); EP(3); EP(4)]/(1+EP(1));
end
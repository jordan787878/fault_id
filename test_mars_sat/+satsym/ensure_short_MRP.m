function [X_short_MRP,aux] = ensure_short_MRP(X)

X_short_MRP = X;

% ensure norm-1 of MRP
MRP = X(1:3);
if(norm(MRP) > 1)
    MRP = -MRP/(norm(MRP)^2);
    aux.MRP_switch = true;
else
    aux.MRP_switch = false;
end
X_short_MRP(1:3) = MRP;

end
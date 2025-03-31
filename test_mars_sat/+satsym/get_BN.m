function BN = get_BN(MRP)
    % BN = EP_to_DCM(MRP_to_EP(MRP));
    den = (1+norm(MRP)^2)^2;
    BN = eye(3) + ( 8*satsym.tilde(MRP)*satsym.tilde(MRP) -4*(1-norm(MRP)^2)*satsym.tilde(MRP) )/ den;
end
% 
% function EP = MRP_to_EP(sigma)
% sigma_square = dot(sigma,sigma);
% EP = [(1-sigma_square)/(1+sigma_square);
%        2*sigma(1)/(1+sigma_square);
%        2*sigma(2)/(1+sigma_square);
%        2*sigma(3)/(1+sigma_square)];
% end
% 
% function DCM = EP_to_DCM(EP)
% b0 = EP(1); 
% b1 = EP(2); 
% b2 = EP(3); 
% b3 = EP(4);
% DCM = [b0^2+b1^2-b2^2-b3^2, 2*(b1*b2+b0*b3), 2*(b1*b3-b0*b2);
%        2*(b1*b2-b0*b3), b0^2-b1^2+b2^2-b3^2, 2*(b2*b3+b0*b1);
%       2*(b1*b3+b0*b2), 2*(b2*b3-b0*b1), b0^2-b1^2-b2^2+b3^2];
% end
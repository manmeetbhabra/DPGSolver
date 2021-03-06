function [Chi_v] = basis_TET_P1(rst)
% rst = rst;

Nb = 4;
Nn = size(rst,1);

a = zeros(Nn,1);
b = zeros(Nn,1);
c = zeros(Nn,1);

for node = 1:Nn
    r = rst(node,1);
    s = rst(node,2);
    t = rst(node,3);    
    
    if (norm(sqrt(6)*t-3) > 1e-12)
        if (norm(2*sqrt(3)*s+sqrt(6)*t-3) > 1e-12)
            a(node,1) = -6*r./(2*sqrt(3)*s+sqrt(6)*t-3);
        end
        b(node,1) = 1/3*(8*sqrt(3)*s./(3-sqrt(6)*t)-1);
    end
    c(node,1) = 1/2*(sqrt(6)*t-1);
end

Chi_v = (4/2^0.25)*...
    [sqrt(1/2)*...
     sqrt(1/2)*...
     sqrt(3/8)+0*a ...
     ...
     sqrt(1/2)*...
     sqrt(1/2)*...
     sqrt(5/8)*1/2*(4*c+2) ...
     ...
     (1-c)*...
     sqrt(1/2)*...
     sqrt(2/2)*1/2.*(3*b+1)*...
     sqrt(5/32) ...
     ...
     (1-b).*(1-c)*...
     sqrt(3/2).*a*...
     sqrt(1/4)*...
     sqrt(5/32)];
 
return



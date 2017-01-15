function y = apply_Bfull(obj,x0)
% The B referred to is the one in Tfull(k) = A - k^2*B.

% split x0 into two pieces
x1 = x0(obj.sc.I1); 
x2 = x0(obj.sc.I2);

x1hat = reshape(x1,obj.dtns.Nr,obj.dtns.Nt);
x1hat(obj.dtns.valRows,:) = 0;
x1hat(obj.dtns.derRows,:) = 0;
x1hat(end,:) = 0;

y = [x1hat(:); x2];
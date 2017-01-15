function x = apply_A_star(obj,x0)
% Does x = dtn.A'*x0

xhat = reshape(x0, obj.Nr, []); % for applying kronecker products

y = obj.localPlacement'*xhat; y = y(:);

x = obj.apply_Aorig_star(x0,obj.VvaluesVec) + ...
    obj.ArowChange'*y;

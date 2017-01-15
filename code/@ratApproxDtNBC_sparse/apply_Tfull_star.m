function y = apply_Tfull_star(obj,x,k)
% does y = rat.Tfull(k)'*x

E = k^2;

% split into two pieces
x1 = x(1:obj.dtns.Nt*obj.dtns.Nr);
x2 = x(obj.dtns.Nt*obj.dtns.Nr+1:end);

% split Ax = [Ax1; Ax2]
Ax1 = obj.dtns.apply_A_star(x1) + obj.apply_A21_star(x2);
Ax2 = obj.apply_A12_star(x1) + conj(obj.A22diag).*x2;

y = [Ax1; Ax2] - conj(E)*obj.apply_Bfull(x);

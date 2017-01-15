function x = apply_Tfull_star_Tfull(obj,x0,k)
% does x = rat.Tfull(k)'*rat.Tfull(k)*x0

x = obj.apply_Tfull_star(obj.apply_Tfull(x0,k),k);

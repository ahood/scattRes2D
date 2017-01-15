function x = apply_B(obj,x0)
% The B referred to is the one in T(k) = A - k^2*B.
% B is the identity with interface condition and boundary condition
% rows zeroed.

x = x0; % B acts like identity
x(obj.globalICBCrows) = 0; % except zeros out certain entries

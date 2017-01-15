function x = apply_Bfull(obj,x0)
% The Bfull is the B in Tfull(k) = A - k^2*B.

% zero boundary and interface condition entries
x = x0;
x(obj.globalICBCrows) = 0;

function y = apply_Uinv(obj,x)
% y = U\x where len(x) = len(obj.BCrows) and x is a vector or
% a matrix considered as a bunch of column vectors.

    y = fft( 1/obj.Nt * diag(exp(-1i*obj.theta.'*(-obj.Nt/2 + 1)))*x );

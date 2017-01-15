function y = apply_U(obj,x)
% y = U*x where len(x) = len(obj.BCrows) and x is a vector or
% a matrix considered as a bunch of column vectors.

    % works, diagonal matrix not formed
    y = obj.Nt*ifft(x);
    d = exp(1i*obj.theta.'*(-obj.Nt/2 + 1));
    for col = 1:size(x,2)
        y(:,col) = d.*y(:,col);
    end

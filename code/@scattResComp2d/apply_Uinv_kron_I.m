function y = apply_Uinv_kron_I(obj,x)
% y = kron(Uinv,I)*x, x = xhat(:)
% If x is a vector, user can pass xhat instead.

    if size(x) == [obj.Nr, obj.Nt] % passed xhat for vector x
        xhat = x;
        yhat = obj.apply_Uinv(xhat.').';
        y = yhat(:);
    else
        y = x; % same shape
        for j = 1:size(x,2)
            xhat = reshape(x(:,j),obj.Nr,obj.Nt);
            yhat = obj.apply_Uinv(xhat.').';
            y(:,j) = yhat(:);
        end
    end
    
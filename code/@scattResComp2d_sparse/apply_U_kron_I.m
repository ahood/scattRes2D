function y = apply_U_kron_I(obj,x)
% y = kron(U,I)*x, x = xhat(:)

    xhat = reshape(x,obj.Nr,obj.Nt);
    yhat = obj.apply_U(xhat.').';
    y = yhat(:);

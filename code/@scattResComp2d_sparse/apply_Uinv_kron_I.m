function y = apply_Uinv_kron_I(obj,x)
% y = kron(Uinv,I)*x, x = xhat(:)

    xhat = reshape(x,obj.Nr,obj.Nt);
    yhat = obj.apply_Uinv(xhat.').';
    y = yhat(:);

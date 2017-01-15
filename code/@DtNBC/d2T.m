function d2Tk = d2T(obj,k)
% second derivative with respect to k
    d = diag( obj.d2DtNcoeffs(obj.Ns,k,obj.r(end)) );
    X = obj.U*d*obj.Uinv;
    Y = zeros(obj.Nr); Y(end,end) = 1;
    d2Ck = -kron(X,Y);

    d2Tk = -2*obj.B + d2Ck;
end
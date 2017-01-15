function Ck = C(obj,k)
% returns the nonlinear part of T
    d = sparse(diag( obj.DtNcoeffs(obj.Ns,k,obj.r(end)) ));
    X = sparse(obj.U*d*obj.Uinv);
    Y = sparse(obj.Nr,obj.Nr); Y(end,end) = 1;
    Ck = -kron(X,Y);
end
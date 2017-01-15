function dTk = dT(obj,k)
% derivative of T w.r.t. k (used for refinement)
    d = diag( obj.dDtNcoeffs(obj.Ns,k,obj.r(end)) );
    X = obj.U*d*obj.Uinv;
    Y = zeros(obj.Nr); Y(end,end) = 1;
    dCk = -kron(X,Y);

    dTk = -2*k*obj.B + dCk;
end
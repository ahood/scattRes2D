function Ck = C_fc(obj,k)
% returns the nonlinear part of T_fc
    d = obj.DtNcoeffs(obj.Ns,k,obj.r(end));
    n = obj.Nr*obj.Nt;
    Ck = sparse(n,n);
    Ck(obj.BCrows,obj.BCrows) = -diag(d);

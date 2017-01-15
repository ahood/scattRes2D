function df = dDtNcoeffs(n,k,R)
% returns the derivative of each DtN coeff w.r.t. k
    fn   = DtNBC.DtNcoeffs(n  ,k,R);
    fnp1 = DtNBC.DtNcoeffs(n+1,k,R);
    fnm1 = DtNBC.DtNcoeffs(n-1,k,R);
    Hn   = besselh(n  ,1,k*R);
    Hnp1 = besselh(n+1,1,k*R);
    Hnm1 = besselh(n-1,1,k*R);
    df = fn.*(1-R*fn)./k + ...
            R*(Hnm1.*fnm1 - Hnp1.*fnp1)/2./Hn;
end
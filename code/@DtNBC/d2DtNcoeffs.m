function d2f = d2DtNcoeffs(n,k,R)
% second deriv with respect to k
    fn   = DtNBC.DtNcoeffs(n  ,k,R);
    fnp1 = DtNBC.DtNcoeffs(n+1,k,R);
    fnm1 = DtNBC.DtNcoeffs(n-1,k,R);
    dfn   = DtNBC.dDtNcoeffs(n  ,k,R);
    dfnp1 = DtNBC.dDtNcoeffs(n+1,k,R);
    dfnm1 = DtNBC.dDtNcoeffs(n-1,k,R);
    Hn   = besselh(n  ,1,k*R);
    Hnp1 = besselh(n+1,1,k*R);
    Hnm1 = besselh(n-1,1,k*R);
    Hnp2 = besselh(n+2,1,k*R);
    Hnm2 = besselh(n-2,1,k*R);
    dHn   = (Hnm1 - Hnp1)*R/2;
    dHnm1 = (Hnm2 - Hn  )*R/2;
    dHnp1 = (Hn   - Hnp2)*R/2;

    dx1 = (dfn./k - fn./k.^2).*(1 - R*fn) - R*fn.*dfn./k;
    x2numer = Hnm1.*fnm1;
    dx2numer = dHnm1.*fnm1 + Hnm1.*dfnm1;
    dx2 = (dx2numer.*Hn - x2numer.*dHn)./Hn.^2;
    x3numer = Hnp1.*fnp1;
    dx3numer = dHnp1.*fnp1 + Hnp1.*dfnp1;
    dx3 = (dx3numer.*Hn - x3numer.*dHn)./Hn.^2;

    d2f = dx1 + (dx2 - dx3)*R/2;
end
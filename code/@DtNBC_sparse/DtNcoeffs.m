function f = DtNcoeffs(n,k,R)
% returns the DtN coefficients f_n(k)
% Asymptotics for large n derived as in my thesis, 
% asymptotics for large k derived from
% http://dlmf.nist.gov/10.17.

    Hn = besselh(n,1,k*R);
    dHn = (besselh(n-1,1,k*R) - besselh(n+1,1,k*R))/2;
    f = k.*dHn./Hn;
    
    if any(isnan(f))
        % bookkeeping for dealing with nans
        if length(n) < length(k)
            n = n*ones(size(k));
        else
            k = k*ones(size(n));
        end
        % nan because n too big 
        idx1 = isnan(f) & (abs(n) > abs(k));
        f(idx1) = -abs(n(idx1))/R + k(idx1).^2*R/4./(abs(n(idx1))-1);
        % nan because k too big
        idx2 = isnan(f) & (abs(k) > abs(n));
        f(idx2) = 1i*k(idx2) - 1/2/R;
    end
end
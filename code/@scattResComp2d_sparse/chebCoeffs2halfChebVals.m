function v = chebCoeffs2halfChebVals(c,parity)
% See halfChebVals2chebCoeffs. This does opposite.

    N = length(c);
    ts = (-N+1:0)*pi/(2*N-1);
    v = cos(ts.'*(0:N-1))*c;

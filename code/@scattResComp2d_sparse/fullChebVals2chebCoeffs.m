function c = fullChebVals2chebCoeffs(v)
% Takes values of an arbitrary function on chebyshev mesh of [-1,1] and
% returns as many chebyshev coefficients.

    N = size(v,1)-1; % N + 1 values
    V = [v; flipud(v(2:end-1,:))]; % extend to vector of 2N values

    d1 = exp( (pi*1i/N)*(-N:  N-1)*(-N+1) ); D1inv = spdiags(1./d1.',0,2*N,2*N);
    d2 = exp( (pi*1i/N)*( 0:2*N-1)*(-N  ) ); D2inv = spdiags(1./d2.',0,2*N,2*N);

    C = D2inv*fft( (D1inv*V)/2/N ); % 2N Fourier coeffs C_{-N+1}, ... , C_0, ..., C_N
    c = C(N:end,:); % take C_0, ..., C_N (N+1 of them)
    c(2:end,:) = 2*c(2:end,:); % scale to get Chebyshev coeffs

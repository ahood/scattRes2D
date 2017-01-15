function v = chebCoeffs2fullChebVals(c)
% Takes Chebyshev coefficients and returns function values at Chebyshev
% grid on [-1,1].

    N = size(c,1)-1; % N + 1 cheb coeffs c_0, ..., c_N
    c(2:end,:) = c(2:end,:)/2; % scale to Fourier coeffs C_0, ..., C_N
    C = [flipud(c(2:end-1,:)); c]; % 2N Fourier coeffs C_{-N+1}, ..., C_0, ..., C_N

    d1 = exp( (pi*1i/N)*(-N:  N-1)*(-N+1) ); D1 = spdiags(d1.',0,2*N,2*N);
    d2 = exp( (pi*1i/N)*( 0:2*N-1)*(-N  ) ); D2 = spdiags(d2.',0,2*N,2*N);

    V = 2*N*D1*ifft(D2*C); % vector of 2N values, first N+1 values on [-1,1]
    v = V(1:N+1,:);
    
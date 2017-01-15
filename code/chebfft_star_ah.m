% CHEBFFT  Chebyshev differentiation via FFT. Simple, not optimal.  
%          If v is complex, delete "real" commands.

function w = chebfft_star_ah(v,a,b) % AH
  N = size(v,1)-1; if N==0, w=0; return, end
  x = cos((0:N)'*pi/N);
  ii = 0:N-1;

  dW = 1i*[ii 0 1-N:-1];
  d1 = ii.^2;
  d2 = 1./sqrt(1-x(2:N).^2);
  d3 = (-1).^(ii+1).*ii.^2;

  % apply first column of Mstar to top row of v
  V1 = (1/N)*[d1'; 0*d1']*v(1,:);
  V1(N+1,:) = V1(N+1,:) + (N/2)*v(1,:);

  % apply middle columns of Mstar to middle rows of v
  tmp = diag(d2')*v(2:N,:);
  tmp2 = [0*tmp(1,:); tmp; 0*v(1:N,:)];
  V2 = -(1/2/N)*diag(dW')*fft(tmp2);
  
  % apply last column of Mstar to last row of v
  V3 = (1/N)*[d3'; 0*d3']*v(N+1,:);
  V3(N+1,:) = V3(N+1,:) + (N/2)*(-1)^(N+1)*v(N+1,:);

  % Mstar*v
  V = V1 + V2 + V3;

  % apply fft_star
  U = 2*N*ifft(V);

  % flip and sum
  U1 = U(1:N+1,:);
  U2 = U(N+2:end,:);
  U2flipped = flipud(U2);
  W = U1;  
  W(2:N,:) = W(2:N,:) + U2flipped;

  % do the scaling
  w = -W;
  w = w*2/(b-a);


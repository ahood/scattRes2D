% CHEBFFT  Chebyshev differentiation via FFT. Simple, not optimal.  
%          If v is complex, delete "real" commands.

function w = chebfft_ah(v,a,b) % AH
  N = size(v,1)-1; if N==0, w=0; return, end
  x = cos((0:N)'*pi/N);
  ii = 0:N-1;
  V = [v; flipud(v(2:N,:))]; % AH - assume v is matrix whose columns are f(x) 
  U = fft(V);
  W = ifft(diag(1i*[ii 0 1-N:-1])*U); % AH - assume v is matrix whose columns are f(x) 
  w = zeros(size(v)); % AH - assume v is matrix whose columns are f(x) 
  w(1,:) = sum(diag(ii.^2)*U(ii+1,:))/N + 0.5*N*U(N+1,:); % AH - assume v is matrix whose columns are f(x) 
  w(2:N,:) = -diag(1./sqrt(1-x(2:N).^2))*W(2:N,:); % AH - assume v is matrix whose columns are f(x) 
  w(N+1,:) = sum(diag((-1).^(ii+1).*ii.^2)*U(ii+1,:))/N + ...
      0.5*(-1)^(N+1)*N*U(N+1,:); % AH - assume v is matrix whose columns are f(x) 
          
  w = -w; % AH - since I order cheb nodes from left to right
  w = w*2/(b-a); % AH - assuming v defined on [a,b], not nec. [-1,1]

% CHEBFFT  Chebyshev differentiation via FFT. Simple, not optimal.  
%          If v is complex, delete "real" commands.
% Modified from chebfft.m in Spectral Methods in Matlab--
% https://people.maths.ox.ac.uk/trefethen/chebfft.m

%   function w = chebfft(v)
function w = chebfft(v,a,b) % AH
%   N = length(v)-1; if N==0, w=0; return, end
  N = size(v,1)-1; if N==0, w=0; return, end
  x = cos((0:N)'*pi/N);
  ii = 0:N-1;
%   v = v(:); % AH - assume v is matrix whose columns are f(x) 
%   V = [v; flipud(v(2:N))];      % transform x -> theta  
  V = [v; flipud(v(2:N,:))]; % AH - assume v is matrix whose columns are f(x) 
%   U = real(fft(V));
%   W = real(ifft(1i*[ii 0 1-N:-1]'.*U));
  U = fft(V);
%   W = ifft(1i*[ii 0 1-N:-1]'.*U);
  W = ifft(diag(1i*[ii 0 1-N:-1])*U); % AH - assume v is matrix whose columns are f(x) 
%   w = zeros(N+1,1);
  w = zeros(size(v)); % AH - assume v is matrix whose columns are f(x) 
%   w(2:N) = -W(2:N)./sqrt(1-x(2:N).^2);    % transform theta -> x     
  w(2:N,:) = -diag(1./sqrt(1-x(2:N).^2))*W(2:N,:); % AH - assume v is matrix whose columns are f(x) 
%   w(1) = sum(ii'.^2.*U(ii+1))/N + .5*N*U(N+1);     
  w(1,:) = sum(diag(ii.^2)*U(ii+1,:))/N + 0.5*N*U(N+1,:); % AH - assume v is matrix whose columns are f(x) 
%   w(N+1) = sum((-1).^(ii+1)'.*ii'.^2.*U(ii+1))/N + ...
%               .5*(-1)^(N+1)*N*U(N+1);
  w(N+1,:) = sum(diag((-1).^(ii+1).*ii.^2)*U(ii+1,:))/N + ...
      0.5*(-1)^(N+1)*N*U(N+1,:); % AH - assume v is matrix whose columns are f(x) 
          
  w = -w; % AH - since I order cheb nodes from left to right
  w = w*2/(b-a); % AH - assuming v defined on [a,b], not nec. [-1,1]

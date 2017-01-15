% CHEB  compute D = differentiation matrix, x = Chebyshev grid
% Modified from cheb.m in Spectral Methods in Matlab--
% https://people.maths.ox.ac.uk/trefethen/cheb.m

  function [D,x] = cheb(N,a,b)
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries
  
  % reorder the nodes from left to right
  P = flipud(eye(N+1));
  D = P*D*P;
  x = P*x;
  
  if nargin == 3
    x = (1-x)*a/2 + (1+x)*b/2;
    D = D*2/(b-a);
  end


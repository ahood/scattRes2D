function truescattn = compute_truescattn(n,k,Vs,s)
% Evaluates n-th fourier coefficient of analytically-computed scattering 
% solution on radial direction mesh associated to scattResComp2d object s, 
% where the associated potential is axisymmetric and piecewise constant in 
% the r direction and the incident wave is e^{ikx}.
% Rs is the list of boundaries (not including zero) and Vs is the values
% of the potential on each subregion.
% BEWARE OF BRANCH CUT.

h = 1e-6*rand + 1i*1e-6*rand;

% list of adjusted frequencies
ks = sqrt(k^2 - Vs);

% n-th fourier coefficient of incident wave and its derivative
inc = @(r) 1i^n*s.J(n,k,r); dinc = @(r) 1i^n*s.dJ(n,k,r);

% set up A and b such that A\b is the coefficients of the true solution
N = length(Vs); A = zeros(2*N); b = zeros(2*N,1);
for j = 2:N % for each interior region besides first
    lt = s.Rs(j-1); rt = s.Rs(j); % left and right endpoints of region
    J = 2*(j-1):2*(j-1)+1; % columns in which region j stuff appears
    Ilt = 2*(j-1)-1:2*(j-1); % rows for interface condition for left  endpoint
    Irt = 2*j    -1:2*j;     %                                  right 
    kj = ks(j); % frequency of waves on this region

    A(Ilt,J) = [  -s.J(n,kj,lt),  -s.Y(n,kj,lt) ; ...
                 -s.dJ(n,kj,lt), -s.dY(n,kj,lt) ];
    A(Irt,J) = [   s.J(n,kj,rt),   s.Y(n,kj,rt) ; ...
                  s.dJ(n,kj,rt),  s.dY(n,kj,rt) ];
end
% region 1
kj = ks(1); rt = s.Rs(1);
A(1:2,1) = [  s.J(n,kj,rt) ; ...
             s.dJ(n,kj,rt) ];
% boundary condition enforcing total wave = incident wave + first kind
% hankel
R = s.Rs(end);
A(end-1:end,end) = [  -s.H(n,k,R) ; ...
                     -s.dH(n,k,R) ];
b(end-1:end) = [  inc(R) ; ...
                 dinc(R) ];
             
% solve
x = A\b;
% x = double(sym(A)\sym(b));
    
% for each region, evaluate piecewise solution for total wave,
%                  subtract incident wave evaluated on that region, 
%                  and append resulting evaluation of scattered wave
truescattn = zeros(sum(s.Nrs),1); 
% first region
I1 = 1:s.Nrs(1); r1 = s.r(I1); k1 = ks(1);
truescattn(I1) = x(1)*s.J(n,k1,r1) - inc(r1);
% all others
for j = 2:N
    Ij = sum(s.Nrs(1:j-1)) + (1:s.Nrs(j)); % indices corresponding to mesh points for j-th region
    rj = s.r(Ij); % mesh on j-th region
    kj = ks(j);
    a = x(2*(j-1)); b = x(2*(j-1)+1); % coeffs of J and Y
    truescattn(Ij) = a*s.J(n,kj,rj) + b*s.Y(n,kj,rj) - inc(rj);
end

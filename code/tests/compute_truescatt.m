function truescatt = compute_truescatt(k,Vs,s,maxn)
% Evaluates true solution to scattering problem on mesh of disc.
% Assumes incident wave is e^{ikx}. BEWARE OF BRANCH CUT.
% Inputs:
%   k  - frequency parameter
%   Vs - values of piecewise constant potential on annular regions given by
%        s.Rs
%   s  - scattResComp2d instance where disc and mesh are specified
%   maxn - fourier coefficients -maxn+1:maxn will be used for computing the
%          true solution.

truescattFourier = zeros(s.Nr,2*maxn);
for n = -maxn+1:maxn; % Nt is even, so number of terms should be too
    truescattFourier(:,n+maxn) = compute_truescattn(n,k,Vs,s);
end
truescatt = s.valuesVecFromFourierVec(truescattFourier(:));

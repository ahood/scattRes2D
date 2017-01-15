function [azns,bzns] = hankel_zeros_small(n,smax,tol,itermax)
% This function computes zeros of the first kind Hankel
% function of order n using the Olver expansion from
% Boro Doering's paper, assuming that the zero is "small".
% 
% Inputs:
%  - n:       order of Hankel function
%  - smax:    the most negative Airy function zeros used in 
%             computing zeros of H_n^1 are a_{smax} and 
%             b_{smax} (see hankel_zeros.pdf notes)
%  - tol:     absolute error tolerance on H_n^1(z_s) in order
%             for z_s to be returned
%  - itermax: maximum number of Newton steps to perform
%
% Outputs:
%  - zns: a vector of numbers such that H_n^1(zns(j)) < tol, all j


if n == 0, azns = []; bzns = []; return; end
if isempty(tol),     tol = 1e-8;   end
if isempty(itermax), itermax = 50; end

sigmatol = tol^2; % just guessing, but probably want this to be a lot more accurate

[as,bs] = airy_zeros(smax);

% first handle as stuff
gs = as*exp(-2*pi*1i/3);
ts = get_ts(gs,n,sigmatol,itermax);
Ws = n*ts;
azns = refine(Ws,n,tol,itermax);

% now to bs
gs = bs*exp(2*pi*1i/3);
ts = get_ts(gs,n,sigmatol,itermax);
Ws = -n*ts;
bzns = refine(Ws,n,tol,itermax);

end

function zns = refine(Ws,n,tol,itermax)
    H = @(z) besselh(n,1,z);
    dH = @(z) (besselh(n-1,1,z)-besselh(n+1,1,z))/2;
    
    resid = H(Ws);
    badinds = abs(resid) > tol;
    iternum = 0;
    while ~isempty(badinds) && iternum < itermax
        Ws = Ws - H(Ws)./dH(Ws);
        
        resid = H(Ws);
        badinds = abs(resid) > tol;
        iternum = iternum + 1;
    end
    zns = Ws( ~badinds );
end

function ts = get_ts(gs,n,sigmatol,itermax)
    rhos = (2/3)*(1/n)*gs.^(3/2);
    sgn = 2*((imag(gs) > 0) - 1/2);
    sigmas = 1.1 + sgn*pi*1i/4; % initial guesses
    sigmas = get_sigma(sigmas,rhos,sigmatol,itermax); % update
    ts = 1./cosh(sigmas);
end

function sigma = get_sigma(sigma,rho,tol,itermax)
    resid = sigma - tanh(sigma) - rho;
    badinds = abs(resid) > tol;
    iternum = 0;
    while ~isempty(badinds) && iternum < itermax
        s = sigma(badinds); r = rho(badinds);
        Q = (s - tanh(s) - r)./tanh(s).^2;
        innerQcoeff = 1./tanh(s) - tanh(s);
        outerQcoeff = 1 + innerQcoeff.*Q;
        s = s - outerQcoeff.*Q;
        sigma(badinds) = s;
        
        resid = sigma - tanh(sigma) - rho;
        badinds = abs(resid) > tol;
        iternum = iternum + 1;
    end
end
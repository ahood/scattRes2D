function zns = hankel_zeros_large(n,smax,tol,itermax)
% This function computes zeros of the first kind Hankel
% function of order n using the McMahon expansion from
% Boro Doering's paper, assuming that the zero is "large".
% 
% Inputs:
%  - n:       order of Hankel function
%  - smax:    the right-most zero we compute is z_{smax}
%  - tol:     absolute error tolerance on H_n^1(z_s) in order
%             for z_s to be returned
%  - itermax: maximum number of Newton steps to perform
%
% Outputs:
%  - zns: a vector of numbers such that H_n^1(zns(j)) < tol, all j


if isempty(tol),     tol     = 1e-8; end
if isempty(itermax), itermax = 50;   end

zns = 0*(1:smax);

% for Newton refinement
H = @(z) besselh(n,1,z);
dH = @(z) (besselh(n-1,1,z)-besselh(n+1,1,z))/2;

% set up D coefficients (they depend on n but not s)
numax = 5;
D = get_D(numax,n);

alpha = 1i*log(2)/2; % figure 3
e = 2*(0:numax-1)+1;
for s = 1:smax
    beta = (s + n/2 - 1/4)*pi - alpha;
    zs = beta - sum( D./beta.^e );
    zs = -conj(zs);
    zs = refine(zs,H,dH,tol,itermax);
    zns(s) = zs;
end

zns = zns( abs(H(zns)) < tol );

end

function zs = refine(zs,H,dH,tol,itermax)
    resid = H(zs);
    badinds = abs(resid) > tol;
    iternum = 0;
    while ~isempty(badinds) && iternum < itermax
        zs = zs - H(zs)/dH(zs);
        
        resid = H(zs);
        badinds = abs(resid) > tol;
        iternum = iternum + 1;
    end
end

function D = get_D(numax,n)
    % initialize
    t0 = 1; t = []; % t shouldn't include t0
    A0 = 1; A = A0; % A should    include A0
    Q = [];         % Q shouldn't include Q0
    C0 = 1; C = C0; % C should    include C0
    D = [];         % D shouldn't include D0
    F = zeros(numax);

    nu = 1;
    t = [t, get_tnu(nu,t0,n)];
    A = [A, get_Anu(nu,A,t)];
    F =    update_F(nu,F,Q,C);
    C = [C, get_Cnu(nu,A(end))];
    D = [D, get_Dnu(nu,D,F,C)];

    for nu = 2:numax
        t = [t, get_tnu(nu,t(end),n)];
        A = [A, get_Anu(nu,A,t)];
        Q = [Q, get_Qnu(   C)];
        F =    update_F(nu,F,Q,C);
        C = [C, get_Cnu(nu,A(end))];
        D = [D, get_Dnu(nu,D,F,C)];
    end
end

function tnu = get_tnu(nu,tprev,n)
    x = nu-1/2;
    tnu = -x*(n^2 - x^2)*tprev/nu;
end

function Anu = get_Anu(nu,A,t)
    l = 0:(nu-1);
    Anu = -sum( (-1).^(nu-l).*A.*flip(t) );  % t should not include t0, A should include A0
end

function Qnu = get_Qnu(C)
    Qnu = sum( C.*flip(C) ); % C0 needs to be included in C
end

function F = update_F(nu,F,Q,C)
    mu = 1;
    l = 1:nu-mu;
    F(nu,mu) = C(end) - Q(l)*F(nu-l,mu);
    for mu = 2:nu
        l = 1:nu-mu;
        F(nu,mu) = F(nu-1,mu-1) - Q(l)*F(nu-l,mu);
    end
end

function Cnu = get_Cnu(nu,Anu)
    Cnu = -Anu/(2*nu-1);
end

function Dnu = get_Dnu(nu,D,F,C)
    l = 1:nu-1;
    Dnu = C(end);
    if ~isempty(l)
        Dnu = Dnu + sum( D(l).*F(nu,l) );
    end
end

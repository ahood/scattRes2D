function k = bordered_newton(k, T, dT, itmax)
% Gets better approximation to eigenvalue of matrix-valued function T(k).
%  Inputs: 
%  - k: initial guess for eigenvalue of T.
%  - T: matrix-valued function T(z).
%  - dT: T'(z).
%  - itmax: number of iterations (default 5).

    if nargin < 4, itmax = 5; end
    
    % setup
    Tk = T(k);
    n = length(Tk);
    a = rand(n,1);
    b = rand(n,1);

    for iter = 1:itmax
        vgk = [Tk, b; a', 0]\[0*b; 1];
        v = vgk(1:end-1); gk = vgk(end);
        dvgk = -[Tk, b; a', 0]\[dT(k)*v; 0];
        dgk = dvgk(end);
        k = k - gk/dgk;
        Tk = T(k);
    end

% Discussion:
% Idea is in the system
% T(k)*v + g(k)*b = 0 (T(k)*v = 0 if and only if g(k) = 0)
% a'*v = 1            (normalization).
%
% In matrix form, [T(k), b; a', 0]*[v; g(k)] = [0, 1].
%
% Given initial guess k, Newton update is k <- k - g(k)/g'(k).
% So we need g(k), g'(k), and incidentally also v.
%
% First of all,
%
% [v; g(k)] = [T(k), b; a', 0]\[0,1].
%
% Derivative of original system wrt k gives
% [T'(k), 0; 0, 0]*[v; g(k)] + [T(k), b; a', 0]*[0; g'(k)] = 0, giving
%
% [0; g'(k)] = [T(k), b; a', 0]\[T'(k)*v; 0].

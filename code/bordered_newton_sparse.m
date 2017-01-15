function k = bordered_newton(k, apply_T, apply_M, apply_dT, apply_dM, n, itermax)
% Gets better approximation to eigenvalue of matrix-valued function T(k).
%  Inputs: 
%  - k: initial guess.
%  - apply_T: apply_T(x,k) performs T(k)*x.
%  - apply_M: apply_M(x,k) performs M(k)*x where M(k) is preconditioner.
%  - apply_dT: apply_dT(x,k) performs T'(k)*x.
%  - apply_dM: apply_dM(x,k) performs dM(k)*x where dM(k) precond for dT.
%  - n: dimension of T(k).

    if nargin < 7: itermax = 5; end

    a = rand(n,1);
    b = rand(n,1);

    function y = apply_Tinv(x,k)
        y = gmres(@(v) apply_T(v,k),x,[],1e-12,n, @(v) apply_M(v,k));
    end
        
    function y = apply_invTba0(x,k)
        % y = [T(k), b; a', 0]\x
        x1 = x(1:n); x2 = x(end);
        Tinvb  = apply_Tinv(b,k);
        Tinvx1 = apply_Tinv(x1,k);
        aTinvb = a'*Tinvb;
        
        y1 = Tinvx1 - Tinvb*(a'*Tinvx1)/aTinvb + Tinvb*x2/aTinvb;
        y2 = a'*Tinvx1/aTinvb - x2/aTinvb;
        y = [y1; y2];
    end

    for iter = 1:itermax
        vgk = apply_invTba0([0*b; 1], k);
        v = vgk(1:end-1); gk = vgk(end);
        zerodgk = apply_invTba0([apply_dT(v,k); 0], k);
        dgk = zerodgk(end);
        k = k - gk/dgk;
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

end
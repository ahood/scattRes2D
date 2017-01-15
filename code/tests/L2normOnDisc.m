function L2norm = L2normOnDisc(f_ontheta,theta)
% Computes the L2-norm, sqrt( int_B(0,R) |f|^2 dA ), of a function f
% whose values along regularly spaced rays are given by chebfuns.
%
% Inputs:
% - theta is a vector of equally spaced points in [0,2*pi],
%   with theta(1) = 0 and theta(end) = 2*pi.
% - f_ontheta{n} is the chebfun representing the value of f on theta(n).
%   We assume f_ontheta{n} is a chebfun on [-R,R] and
%   f_ontheta{1} = f_ontheta{end}.

% Define g(t) = int_ 0^R |f(r,t)|^2  r   dr
%             = int_ 0^R | f(r,t)^2  r | dr
% given the way sum works on chebfuns.

% Compute g(theta(n)) for each n.
R = f_ontheta{1}.domain(end);
rfun = chebfun(@(r) r, [-R,R], 'splitting', 'on');

g_ontheta = 0*theta;
for n = 1:length(theta)
    g_ontheta(n) = sum(abs( f_ontheta{n}.^2 .* rfun ), 0, R);
end

% The L2-norm of f is sqrt( int_0^2*pi g(t) dt ). 
L2norm = sqrt( trapz(theta,g_ontheta) );
function [D2r,Dr,r] = radialCheb(N,R1,R2,symmetry)
% takes regular cheb grid on [-1,1] and diff matrix and enforces that it's
% only acting on symmetric functions, i.e. symmetric vectors f =
% [flipud(f2); f2]. Then df2 = D21*flipud(f2) + D22*f2, or df2 =
% (fliplr(D21) + D22)*f2. Take Dr = fliplr(D21)+D22.

if nargin < 4, symmetry = 1; end

[D,x] = cheb(2*N-1); % full grid on [-1,1]

half = length(x)/2; I1 = 1:half; I2 = half+1:2*half;
r = x(I2); % half grid on [0,1]
D2 = D*D;

if symmetry
    Dr  =  fliplr(D (I2,I1)) + D (I2,I2);
    D2r =  fliplr(D2(I2,I1)) + D2(I2,I2);
else
    Dr  = -fliplr(D (I2,I1)) + D (I2,I2);
    D2r = -fliplr(D2(I2,I1)) + D2(I2,I2);
end

% scale
L = R2-R1;
r = L*r + R1;
Dr = Dr/L;
D2r = D2r/L/L;
function [as,bs] = airy_zeros(smax)

% initial guesses
ainds = 1:smax;
binds = 0:smax;
as = -(1.5*pi*(ainds-0.25)).^(2/3);
bs = -(1.5*pi*(binds+0.25)).^(2/3);

% clean up with Newton refinement
for ii = 1:5
    as = as - airy(0,as)./airy(1,as);
    bs = bs - airy(2,bs)./airy(3,bs);
end

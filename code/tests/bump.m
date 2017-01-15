function bumpx = bump(x)
% bump function supported on [-1,1] with max height 1.

    bumpx = (abs(x) < 1).*exp(-1./(1-x.^2));
    bumpx( isnan(bumpx) ) = 0;


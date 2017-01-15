function c = halfChebVals2chebCoeffs(v,parity)
% Takes the values of a function on a half Chebyshev mesh of [0,1] and
% returns as many Chebyshev coefficients. Function symmetry specified as
% 'even' or 'odd' through the parity variable. The variable v can be a
% vector of values or a matrix whose columns are vectors of values.

    % extend to values on cheb mesh of [-1,1] with 2N points
    if strcmp(parity,'even')
        V = [ flipud(v); v]; % extend     symmetrically to [-1,0]
    elseif strcmp(parity,'odd')
        V = [-flipud(v); v]; % extend antisymmetrically to [-1,0]
    end

    C = scattResComp2d.fullChebVals2chebCoeffs(V); % 2N Chebyshev coeffs
    
    % take the coefficients of the even or odd Cheb polys
    if strcmp(parity,'even')
        c = C(1:2:end,:);
    else
        c = C(2:2:end,:);
    end

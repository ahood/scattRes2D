function fourierVec = fourierVecFromChebVec(obj,chebVec,deriv)
% See chebVecFromFourierVec. This does the opposite.

    % symmetric -> asymmetric under deriv, and vice versa
    if nargin < 3, deriv = 0; end

    chebRect = reshape(chebVec,obj.Nr,obj.Nt);
    fourierRect = chebRect;

    % last fourier coefficient has even index if deriv and Nt/2 have same
    % parity, odd index otherwise. So depends on deriv + Nt/2.
    if mod(obj.Nt/2 + deriv,2) == 0 % halfNt + deriv even
        % last column is even
        fourierRect(:,2:2:end) = obj.chebCoeffs2radialVals(chebRect(:,2:2:end),'even');
        fourierRect(:,1:2:end) = obj.chebCoeffs2radialVals(chebRect(:,1:2:end),'odd' );
    else % halfNt + deriv odd
        % last column is odd
        fourierRect(:,2:2:end) = obj.chebCoeffs2radialVals(chebRect(:,2:2:end),'odd' );
        fourierRect(:,1:2:end) = obj.chebCoeffs2radialVals(chebRect(:,1:2:end),'even');
    end
    
    fourierVec = fourierRect(:);  
    
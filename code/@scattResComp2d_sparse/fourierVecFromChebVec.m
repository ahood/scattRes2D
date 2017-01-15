function fourierVec = fourierVecFromChebVec(obj,chebVec)
% See chebVecFromFourierVec. This does the opposite.

    chebRect = reshape(chebVec,obj.Nr,obj.Nt);
    fourierRect = chebRect;
    
    halfNt = obj.Nt/2;
    if mod(obj.Nt/2,2) == 0 % halfNt even
        % last column is even
        fourierRect(:,2:2:end) = obj.chebCoeffs2radialVals(chebRect(:,2:2:end),'even');
        fourierRect(:,1:2:end) = obj.chebCoeffs2radialVals(chebRect(:,1:2:end),'odd' );
    else % halfNt odd
        % last column is odd
        fourierRect(:,2:2:end) = obj.chebCoeffs2radialVals(chebRect(:,2:2:end),'odd' );
        fourierRect(:,1:2:end) = obj.chebCoeffs2radialVals(chebRect(:,1:2:end),'even');
    end
    
    fourierVec = fourierRect(:);  
    
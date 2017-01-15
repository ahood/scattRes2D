function chebVec = chebVecFromFourierVec(obj,fourierVec)
% Takes a vector of the form 
%    [f_{-maxn+1}(obj.r); f_{-maxn+2}(obj.r); ... ; f_{maxn}(obj.r)]
% (see fourierVecFromValuesVec)
% where obj.r is a piecewise radial mesh, and returns
%    [c_{-maxn+1}; c_{-maxn+2}; ... ; f_{maxn}],
% each c_j corresponding to f_j(obj.r) and consisting of Chebyshev
% coefficients. 
% The radial mesh obj.r consists of the concatenation of Chebyshev meshes
% defined on [0,obj.Rs(1)], [obj.Rs(1),obj.Rs(2)], ... separately.
% Therefore each f_j is piecewise defined, and c_j is the concatenation of
% the Chebyshev coefficients for f_j on [0,obj.Rs(1)],
% [obj.Rs(1),obj.Rs(2)], ... separately.

    fourierRect = reshape(fourierVec,obj.Nr,obj.Nt);
    chebRect = fourierRect;
    
    halfNt = obj.Nt/2;
    if mod(obj.Nt/2,2) == 0 % halfNt even
        % last column is even
        chebRect(:,2:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,2:2:end),'even');
        chebRect(:,1:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,1:2:end),'odd' );
    else % halfNt odd
        % last column is odd
        chebRect(:,2:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,2:2:end),'odd' );
        chebRect(:,1:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,1:2:end),'even');
    end
    
    chebVec = chebRect(:);

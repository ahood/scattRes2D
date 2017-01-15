function chebVec = apply_W(obj,fourierVec,deriv)
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
%
% Also works if input is a matrix whose columns are considered fourierVecs.

    % symmetric -> asymmetric under deriv, and vice versa
    if nargin < 3, deriv = 0; end

    X = fourierVec; % neutral name since could be matrix or vector
    X_reshaped = reshape(X, 2*obj.Nr, []);
    Y_reshaped = X_reshaped; % to get right size
    if mod(obj.Nt/2 + deriv,2) == 0 % last col even
        odd_rows  = 1:obj.Nr;  even_rows = obj.Nr+1:2*obj.Nr;
    else
        even_rows = 1:obj.Nr;  odd_rows  = obj.Nr+1:2*obj.Nr;
    end
    Y_reshaped( odd_rows,:) = obj.radialVals2chebCoeffs(X_reshaped( odd_rows,:), 'odd');
    Y_reshaped(even_rows,:) = obj.radialVals2chebCoeffs(X_reshaped(even_rows,:),'even');
    Y = reshape(Y_reshaped, obj.Nr*obj.Nt, []);
    chebVec = Y; % columns are chebVecs

%     fourierRect = reshape(fourierVec,obj.Nr,obj.Nt);
%     chebRect = fourierRect;
% 
%     % last fourier coefficient has even index if deriv and Nt/2 have same
%     % parity, odd index otherwise. So depends on deriv + Nt/2.    
%     if mod(obj.Nt/2 + deriv,2) == 0 % halfNt + deriv even
%         % last column is even
%         chebRect(:,2:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,2:2:end),'even');
%         chebRect(:,1:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,1:2:end),'odd' );
%     else % halfNt + deriv odd
%         % last column is odd
%         chebRect(:,2:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,2:2:end),'odd' );
%         chebRect(:,1:2:end) = obj.radialVals2chebCoeffs(fourierRect(:,1:2:end),'even');
%     end
%     
%     chebVec = chebRect(:);

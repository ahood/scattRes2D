function chebVec = chebVecFromFourierVec(obj,fourierVec,deriv)
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

    chebVec = obj.apply_W(fourierVec,deriv);

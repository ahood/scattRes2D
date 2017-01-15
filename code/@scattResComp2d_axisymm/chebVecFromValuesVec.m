function chebVec = chebVecFromValuesVec(obj,valuesVec,deriv)

    % symmetric -> asymmetric under deriv, and vice versa
    if nargin < 3, deriv = 0; end

    fourierVec = obj.fourierVecFromValuesVec(valuesVec);
    chebVec = obj.chebVecFromFourierVec(fourierVec,deriv);
    
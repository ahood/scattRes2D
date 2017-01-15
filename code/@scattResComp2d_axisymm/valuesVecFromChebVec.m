function valuesVec = valuesVecFromChebVec(obj,chebVec,deriv)

    % symmetric -> asymmetric under deriv, and vice versa
    if nargin < 3, deriv = 0; end

    fourierVec = obj.fourierVecFromChebVec(chebVec,deriv);
    valuesVec = obj.valuesVecFromFourierVec(fourierVec);
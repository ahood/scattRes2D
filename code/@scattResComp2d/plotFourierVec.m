function u = plotFourierVec(obj,fourierVec,titlestr,part)
% Plots the function represented with fourierVec. See
% valuesVecFromFourierVec.m for what a fourierVec is.

    if nargin < 4, part = @real; end
    
    valuesVec = obj.valuesVecFromFourierVec(fourierVec);
    u = obj.plotValuesVec(valuesVec,titlestr,part);

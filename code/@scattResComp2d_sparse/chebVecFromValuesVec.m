function chebVec = chebVecFromValuesVec(obj,valuesVec)

    fourierVec = obj.fourierVecFromValuesVec(valuesVec);
    chebVec = obj.chebVecFromFourierVec(fourierVec);
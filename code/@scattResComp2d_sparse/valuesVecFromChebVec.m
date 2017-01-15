function valuesVec = valuesVecFromChebVec(obj,chebVec)

    fourierVec = obj.fourierVecFromChebVec(chebVec);
    valuesVec = obj.valuesVecFromFourierVec(fourierVec);
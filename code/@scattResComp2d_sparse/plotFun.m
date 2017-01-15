function u = plotFun(obj,f,titlestr,part,coords)
    if iscell(f)
        valuesVec = obj.valuesVecFromFunCellArray(f,coords);
    else
        valuesVec = obj.valuesVecFromFun(f,coords);
    end
    u = obj.plotValuesVec(valuesVec,titlestr,part);

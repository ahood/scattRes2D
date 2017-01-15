function RHS = RHSfromFun_fc(obj,incfun)
% Converts the values of the RHS in the scattering problem in the "values"
% basis into the "Fourier" basis.

    if isnumeric(incfun)
        k = incfun;
        incfun = @(x,y) exp(1i*k*x);
        incfunVals = obj.valuesVecFromFun(incfun,'rect');
    else
        incfunVals = obj.valuesVecFromFun(incfun,obj.coords);
    end

    Vvals = obj.valuesVecFromFunCellArray(obj.Vs,obj.coords);
    RHS = obj.fourierVecFromValuesVec(-incfunVals.*Vvals);
    
    % zero out rows associated to interface and boundary conditions
    RHS = reshape(RHS,obj.Nr,[]);
    RHS(obj.valRows,:) = 0;
    RHS(obj.derRows,:) = 0;
    RHS(end        ,:) = 0;
    RHS = RHS(:);    

function RHS = RHSfromFun_cc(obj,incfun)
% Converts the values of the RHS in the scattering problem in the "values"
% basis into the "Chebyshev" basis.

    if isnumeric(incfun)
        k = incfun;
        incfun = @(x,y) exp(1i*k*x);
        incfunVals = obj.valuesVecFromFun(incfun,'rect');
    else
        incfunVals = obj.valuesVecFromFun(incfun,obj.coords);
    end

    Vvals = obj.valuesVecFromFunCellArray(obj.Vs,obj.coords);
    RHS = obj.chebVecFromValuesVec(-incfunVals.*Vvals);
    
    % zero out rows associated to interface and boundary conditions
    RHS = reshape(RHS,obj.Nr,[]);
    RHS(obj.valRows_cc,:) = 0;
    RHS(obj.derRows_cc,:) = 0;
    RHS(end        ,:) = 0;
    RHS = RHS(:);    

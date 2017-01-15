function RHS = RHSfromFun_full(obj,incfun)
% Takes incfun (defined in same coordinates as V)
% and sets up vector of values of righthand side 
%    -V*incfun 
% for use in scattering computation, with
% entries in boundary and interface condition rows zeroed.
% If user passes a number for incfun, assumes number passed is
% wavenumber k and uses default incfun = exp(1i*k*x).

    if isnumeric(incfun)
        k = incfun;
        incfun = @(x,y) exp(1i*k*x);
        incfunVals = obj.valuesVecFromFun(incfun,'rect');
    else
        incfunVals = obj.valuesVecFromFun(incfun,obj.coords);
    end

    RHS = -incfunVals.*obj.VvaluesVec;
    
    % zero out rows associated to interface and boundary conditions
    RHS = reshape(RHS,obj.Nr,[]);
    RHS(obj.valRows,:) = 0;
    RHS(obj.derRows,:) = 0;
    RHS(end        ,:) = 0;
    RHS = RHS(:);
end
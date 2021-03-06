function RHS = RHSfromFun(obj,incfun)
% Takes incfun (defined in same coordinates as V)
% and sets up vector of values of righthand side 
%    -V*incfun 
% on mesh of potential support disc (excluding PML region)
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

    Vvals = obj.valuesVecFromFunCellArray(obj.Vs,obj.coords);
    RHS = -incfunVals.*Vvals;
    
    % zero out rows associated to interface and boundary conditions
    RHS = reshape(RHS,obj.Nr,[]);
    RHS(obj.valRows,:) = 0;
    RHS(obj.derRows,:) = 0;
    RHS(end        ,:) = 0;
    
    % get rid of part associated to PML region
    RHS(end-obj.Nrs(end)+1:end,:) = [];
    
    RHS = RHS(:);
end
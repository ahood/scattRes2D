function scattValuesVec = solve(obj,k,incfun)
% Solve the scattering problem to get values of scattered wave on mesh of
% B(0,R).
    if nargin < 3
        RHS = obj.RHSfromFun(k);
    elseif isa(incfun,'numeric')
        RHS = incfun;
    else
        RHS = obj.RHSfromFun(incfun);
    end
    scattValuesVec = obj.T(k)\RHS;

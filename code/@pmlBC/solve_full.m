function scattValuesVec_ext = solve_full(obj,k,incfun)
% solve the scattering problem on disc including PML region

    if nargin < 3
        RHS = obj.RHSfromFun_full(k);
    elseif isvector(incfun)
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_full(incfun);
    end
    scattValuesVec_ext = obj.Tfull(k)\RHS;            

function scattValuesVec_ext = solve_full_cc(obj,k,incfun)
% solve the scattering problem on disc including PML region

    if nargin < 3
        RHS = obj.RHSfromFun_full_cc(k);
    elseif isvector(incfun)
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_full_cc(incfun);
    end
    scattValuesVec_ext = obj.Tfull_cc(k)\RHS;            

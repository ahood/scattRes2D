function scattValuesVec_ext = solve_full_fc(obj,k,incfun)
% solve the scattering problem on disc including PML region

    if nargin < 3
        RHS = obj.RHSfromFun_full_fc(k);
    elseif isvector(incfun)
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_full_fc(incfun);
    end

    scattValuesVec_ext = obj.Tfull_fc(k)\RHS;            

function scattValuesVec_ext = solve_full_cc(obj,k,incfun)
% solve the scattering problem on disc including PML region

    if nargin < 3
        RHS = obj.RHSfromFun_full_cc(k);
    elseif isvector(incfun)
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_full_cc(incfun);
    end
    obj.Tfull_cc(k);
    RHSrect = reshape(RHS,obj.Nr,obj.Nt);
    scattValuesRect = RHSrect;
    for j = 1:obj.Nt
        scattValuesRect(:,j) = obj.Tfull_cc_blks{j}\RHSrect(:,j);
    end
    scattValuesVec_ext = scattValuesRect(:);

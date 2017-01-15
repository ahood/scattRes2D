function scattValuesVec_ext = solve_full_fc(obj,k,incfun)
% solve the scattering problem on disc including PML region

    if nargin < 3
        RHS = obj.RHSfromFun_full_fc(k);
    elseif isvector(incfun)
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_full_fc(incfun);
    end

    obj.Tfull_fc(k);
    RHSrect = reshape(RHS,obj.Nr,obj.Nt);
    scattValuesRect = RHSrect;
    for j = 1:obj.Nt
        scattValuesRect(:,j) = obj.Tfull_fc_blks{j}\RHSrect(:,j);
    end
    scattValuesVec_ext = scattValuesRect(:);

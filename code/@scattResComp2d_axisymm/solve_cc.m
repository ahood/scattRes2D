function scattChebVec = solve_cc(obj,k,incfun)
% Solve the scattering problem to get values of scattered wave on mesh of
% B(0,R).
    if nargin < 3
        RHS = obj.RHSfromFun_cc(k);
    elseif isa(incfun,'numeric')
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_cc(incfun);
    end
    obj.T_cc(k);
    RHSrect = reshape(RHS,[],obj.Nt);
    scattChebRect = RHSrect;
    for j = 1:obj.Nt
        scattChebRect(:,j) = obj.T_cc_blks{j}\RHSrect(:,j);
    end
    scattChebVec = scattChebRect(:);

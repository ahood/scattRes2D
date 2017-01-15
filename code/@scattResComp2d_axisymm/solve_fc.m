function scattFourierVec = solve_fc(obj,k,incfun)
% Solve the scattering problem to get values of scattered wave on mesh of
% B(0,R).
    if nargin < 3
        RHS = obj.RHSfromFun_fc(k);
    elseif isa(incfun,'numeric')
        RHS = incfun;
    else
        RHS = obj.RHSfromFun_fc(incfun);
    end

    RHSrect = reshape(RHS,[],obj.Nt);
    scattFourierRect = RHSrect;
    obj.T_fc(k);
    for j = 1:obj.Nt
        scattFourierRect(:,j) = obj.T_fc_blks{j}\RHSrect(:,j);
    end
    scattFourierVec = scattFourierRect(:);

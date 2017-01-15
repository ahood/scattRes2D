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
    scattFourierVec = obj.T_fc(k)\RHS;
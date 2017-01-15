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
    scattChebVec = obj.T_cc(k)\RHS;

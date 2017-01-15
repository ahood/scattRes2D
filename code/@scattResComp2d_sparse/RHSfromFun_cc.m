% right-hand side for scattering problem (Cheb basis)
function RHS = RHSfromFun_cc(obj,Vs,incfun,coords)
    % Takes incfun and V and sets up vector version of
    % righthand side -V*incfun for use in scattering computation.
    % 1) make fourierVec
    % 2) set rows where BCs are enforced to zero

    % get the values vector and then convert to fourierVec
    incfunVals = obj.valuesVecFromFun(incfun,coords);
    Vvals = obj.valuesVecFromFunCellArray(Vs,coords);
    RHSvals = -incfunVals.*Vvals;
    RHS = obj.fourierVecFromValuesVec(RHSvals);
    
    % convert fourierVec to cheb basis
    if mod(obj.Ns,2)
        W1 = obj.W_as; W2 = obj.W_s;
    else
        W1 = obj.W_s;  W2 = obj.W_as;
    end
        % scale the block rows
    for n = 1:obj.Nt
        In = (n-1)*obj.Nr + 1:n*obj.Nr; % block row indices
        if mod(n,2), RHS(In) = W1\RHS(In); else RHS(In) = W2\RHS(In); end
    end
    
    % zero out rows associated to interface and boundary conditions
    RHS = reshape(RHS,obj.Nr,[]);
    RHS(obj.valRows_cc,:) = 0;
    RHS(obj.derRows_cc,:) = 0;
    RHS(end        ,:) = 0;
    RHS = RHS(:);
end
function enforceDE_fc(obj)
% Expresses A and B in the fourier basis as A_fc and B_fc.

    R = diag(1./obj.r);
    R2 = diag(1./obj.r.^2);
    if isnumeric(obj.Vs)
        VvaluesVec = obj.Vs;
    else
        VvaluesVec = obj.valuesVecFromFunCellArray(obj.Vs,obj.coords);
    end

    % blocks of radial part of Laplacian acting in Fourier space
    Lr_odd  = obj.Drs_even*obj.Drs_odd  + R*obj.Drs_odd;
    Lr_even = obj.Drs_odd *obj.Drs_even + R*obj.Drs_even;
    
    if mod(obj.Nt/2,2) == 0 % last is even
        Lr1 = Lr_odd;  Lr2 = Lr_even;
    else
        Lr1 = Lr_even; Lr2 = Lr_odd;
    end
    
    % set up -Laplacian part of A_fc
    A_fc = zeros(obj.Nr*obj.Nt);
    for j = 1:2:obj.Nt-1
        rows = (j-1)*obj.Nr + 1 : j*obj.Nr; cols = rows;
        A_fc(rows,cols) = -(Lr1 - obj.Ns(j  )^2 * R2);
        rows = j*obj.Nr + 1 : (j+1)*obj.Nr; cols = rows;
        A_fc(rows,cols) = -(Lr2 - obj.Ns(j+1)^2 * R2);
    end

    % add on the potential term
    n = length(VvaluesVec);
    Vdiag = spdiags(VvaluesVec,0,n,n);
    pot_term = obj.apply_U_kron_I_on_right(obj.apply_Uinv_kron_I(Vdiag));
    A_fc = A_fc + pot_term;
        
    obj.A_fc = A_fc;
    obj.B_fc = eye(obj.Nr*obj.Nt);

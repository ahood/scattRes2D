function enforceDE_cc(obj)
% While A and B map values to values, A_cc and B_cc map Chebyshev
% coefficients to Chebyshev coefficients. 

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
        f1 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_odd,  'odd' );
        f2 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_even, 'even');
    else
        Lr1 = Lr_even; Lr2 = Lr_odd;
        f1 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_even, 'even');
        f2 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_odd,  'odd' );
    end
    
    % set up -Laplacian part of A_fc
    A_cc = zeros(obj.Nr*obj.Nt);
    for j = 1:2:obj.Nt-1
        rows = (j-1)*obj.Nr + 1 : j*obj.Nr; cols = rows;
        A_cc(rows,cols) = -f1(Lr1 - obj.Ns(j  )^2 * R2);
        rows = j*obj.Nr + 1 : (j+1)*obj.Nr; cols = rows;
        A_cc(rows,cols) = -f2(Lr2 - obj.Ns(j+1)^2 * R2);
    end
    
    % add on the potential term
    n = length(VvaluesVec);
    Vdiag = spdiags(VvaluesVec,0,n,n); % in values basis
    pot_term = obj.apply_U_kron_I_on_right(obj.apply_Uinv_kron_I(Vdiag)); % convert to fourier basis
    pot_term = obj.apply_W(obj.apply_Winv_on_right(pot_term)); % convert to cheb basis
    A_cc = A_cc + pot_term;
    
    obj.A_cc = A_cc;
    obj.B_cc = eye(obj.Nr*obj.Nt);

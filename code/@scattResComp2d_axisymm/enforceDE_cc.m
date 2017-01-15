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
    Vdiag = diag(VvaluesVec(1:obj.Nr));
    
    % blocks of radial part of Laplacian acting in Fourier space
    Lr_odd  = obj.Drs_even*obj.Drs_odd  + R*obj.Drs_odd;
    Lr_even = obj.Drs_odd *obj.Drs_even + R*obj.Drs_even;
    
    if mod(obj.Nt/2,2) == 0 % last one is even
        Lr1 = Lr_odd;  Lr2 = Lr_even;
        f1 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_odd,  'odd' );
        f2 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_even, 'even');
    else
        Lr1 = Lr_even; Lr2 = Lr_odd;
        f1 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_even, 'even');
        f2 = @(X) obj.radialVals2chebCoeffs(X*obj.Winv_odd,  'odd' );
    end
    
    for j = 1:2:obj.Nt-1
        obj.A_cc_blks{j  } = f1(-(Lr1 - obj.Ns(j  )^2 * R2) + Vdiag);
        obj.A_cc_blks{j+1} = f2(-(Lr2 - obj.Ns(j+1)^2 * R2) + Vdiag);
    end
    for j = 1:obj.Nt
        obj.B_cc_blks{j} = speye(obj.Nr);
    end

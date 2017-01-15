function enforceDE_full_cc(obj)
% While A and B map values to values, A_cc and B_cc map Chebyshev
% coefficients to Chebyshev coefficients. 

    R2 = diag(1./obj.s.^2);
    VvaluesVec = obj.valuesVecFromFunCellArray(obj.Vs,obj.coords);
    Vdiag = diag(VvaluesVec(1:obj.Nr));

    % blocks of radial part of Laplacian acting in Fourier space
    Lr_odd  =  diag(1./obj.s./obj.ds_short - obj.d2s_short./obj.ds_short.^3)*...
               (obj.Dr_I - obj.Dr_J) + ...
               diag(1./obj.ds_short.^2)*(obj.Drr_I - obj.Drr_J);
    Lr_even =  diag(1./obj.s./obj.ds_short - obj.d2s_short./obj.ds_short.^3)*...
               (obj.Dr_I + obj.Dr_J) + ...
               diag(1./obj.ds_short.^2)*(obj.Drr_I + obj.Drr_J);
           
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
        obj.Afull_cc_blks{j  } = -f1(Lr1 - obj.Ns(j  )^2 * R2) + Vdiag;
        obj.Bfull_cc_blks{j  } = speye(obj.Nr);
        obj.Afull_cc_blks{j+1} = -f2(Lr2 - obj.Ns(j+1)^2 * R2) + Vdiag;
        obj.Bfull_cc_blks{j+1} = speye(obj.Nr);
    end

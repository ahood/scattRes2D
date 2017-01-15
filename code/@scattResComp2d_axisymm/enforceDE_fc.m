function enforceDE_fc(obj)
% Expresses A and B in the fourier basis as A_fc and B_fc.

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
    else
        Lr1 = Lr_even; Lr2 = Lr_odd;
    end
    
    for j = 1:2:obj.Nt-1
        obj.A_fc_blks{j  } = -(Lr1 - obj.Ns(j  )^2 * R2) + Vdiag;
        obj.B_fc_blks{j  } = speye(obj.Nr);
        obj.A_fc_blks{j+1} = -(Lr2 - obj.Ns(j+1)^2 * R2) + Vdiag;
        obj.B_fc_blks{j+1} = speye(obj.Nr);
    end

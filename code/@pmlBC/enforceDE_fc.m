function enforceDE_fc(obj)
% Expresses A and B in the fourier basis as A_fc and B_fc.

    R = diag(1./obj.s);
    VvaluesVec = obj.valuesVecFromFunCellArray(obj.Vs,obj.coords);

    % blocks of radial part of Laplacian acting in Fourier space
    Lr_odd  =  diag(1./obj.s./obj.ds_short - obj.d2s_short./obj.ds_short.^3)*...
               (obj.Dr_I - obj.Dr_J) + ...
               diag(1./obj.ds_short.^2)*(obj.Drr_I - obj.Drr_J);
    Lr_even =  diag(1./obj.s./obj.ds_short - obj.d2s_short./obj.ds_short.^3)*...
               (obj.Dr_I + obj.Dr_J) + ...
               diag(1./obj.ds_short.^2)*(obj.Drr_I + obj.Drr_J);

    % set up -Laplacian part of A_fc
    A_fc = zeros(obj.Nr*obj.Nt);
    if mod(obj.Nt/2,2) == 0 % Nt/2 even
        for j = 1:2:obj.Nt-1 % odd blocks
            rows = (j-1)*obj.Nr + 1 : j*obj.Nr; cols = rows;
            A_fc(rows,cols) = -(Lr_odd  - obj.Ns(j)^2 * R^2);
        end
        for j = 2:2:obj.Nt % even blocks
            rows = (j-1)*obj.Nr + 1 : j*obj.Nr; cols = rows;
            A_fc(rows,cols) = -(Lr_even - obj.Ns(j)^2 * R^2);
        end
    else % Nt/2 odd
        for j = 2:2:obj.Nt % odd blocks
            rows = (j-1)*obj.Nr + 1 : j*obj.Nr; cols = rows;
            A_fc(rows,cols) = -(Lr_odd  - obj.Ns(j)^2 * R^2);
        end
        for j = 1:2:obj.Nt-1 % even blocks
            rows = (j-1)*obj.Nr + 1 : j*obj.Nr; cols = rows;
            A_fc(rows,cols) = -(Lr_even - obj.Ns(j)^2 * R^2);
        end
    end

    % add on the potential term
    n = length(VvaluesVec);
    Vdiag = spdiags(VvaluesVec,0,n,n);
    pot_term = obj.apply_U_kron_I_on_right(obj.apply_Uinv_kron_I(Vdiag));
    A_fc = A_fc + pot_term;
        
    obj.A_fc = A_fc;
    obj.B_fc = eye(obj.Nr*obj.Nt);

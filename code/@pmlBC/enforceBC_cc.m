function enforceBC_cc(obj)
% Enforce Dirichlet boundary conditions in Chebyshev basis.

    % zero out the boundary condition rows
    A_cc = reshape(obj.A_cc, obj.Nr, []); A_cc(end,:) = 0; 
    B_cc = reshape(obj.B_cc, obj.Nr, []); B_cc(end,:) = 0;
    A_cc = reshape(A_cc, obj.Nr*obj.Nt, []);
    B_cc = reshape(B_cc, obj.Nr*obj.Nt, []);

    % put last row of Winv blocks in A_cc
    if mod(obj.Nt/2,2) == 0 % last block is even
        Winv1 = obj.Winv_odd;  Winv2 = obj.Winv_even;
    else
        Winv1 = obj.Winv_even; Winv2 = obj.Winv_odd;
    end

    for j = 1:2:obj.Nt
        row = j*obj.Nr;
        cols = row - obj.Nr + 1 : row;
        A_cc(row,cols) = Winv1(end,:);
        row = (j+1)*obj.Nr;
        cols = row - obj.Nr + 1 : row;
        A_cc(row,cols) = Winv2(end,:);
    end

    obj.A_cc = A_cc; obj.B_cc = B_cc;

function enforceBC_cc(obj)
% enforces boundary conditions in cheb basis

    A = obj.A_cc; B = obj.B_cc;
 
    % zero out all rows where boundary conditions will be placed
    A(obj.BCrows,:) = 0; B(obj.BCrows,:) = 0;
    
    % boundary conditions will be enforced on block rows
    for blk = 1:obj.Nt
        I = (blk-1)*obj.Nr + 1:blk*obj.Nr; % indices for this block
        A(I,I) = obj.enforceBlockBC_cc(blk,A(I,I));
    end
    obj.A_cc = A; obj.B_cc = B;

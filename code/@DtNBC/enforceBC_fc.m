function enforceBC_fc(obj)
    A = obj.A_fc; B = obj.B_fc;
 
    % zero out all rows where boundary conditions will be placed
    A(obj.BCrows,:) = 0; B(obj.BCrows,:) = 0;
    
    % boundary conditions will be enforced on block rows
    for blk = 1:obj.Nt
        I = (blk-1)*obj.Nr + 1:blk*obj.Nr; % indices for this block
        A(I,I) = obj.enforceBlockBC_fc(blk,A(I,I));
    end
    obj.A_fc = A; obj.B_fc = B;

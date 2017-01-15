function enforceBC(obj)
    A = obj.A; B = obj.B;
 
    % zero out all rows where boundary conditions will be placed
    A(obj.BCrows,:) = 0; B(obj.BCrows,:) = 0;
    
    % boundary conditions will be enforced on block rows
    for blk = 1:obj.Nt
        I = (blk-1)*obj.Nr + 1:blk*obj.Nr; % indices for this block
        A(I,:) = obj.enforceBlockBC(blk,A(I,:));
    end
    obj.A = A; obj.B = B;


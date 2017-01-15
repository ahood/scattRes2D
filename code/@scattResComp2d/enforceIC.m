function enforceIC(obj)
    A = obj.A; B = obj.B;
    
    % zero out all rows where interface conditions will be placed
    A = reshape(A,obj.Nr,[]); B = reshape(B,obj.Nr,[]); % put all blocks in one row
    A(obj.valRows,:) = 0; B(obj.valRows,:) = 0; % zero out continuity condition rows
    A(obj.derRows,:) = 0; B(obj.derRows,:) = 0; %          differentiability
    A = reshape(A,obj.Nr*obj.Nt,[]); B = reshape(B,obj.Nr*obj.Nt,[]);
    
    % interface conditions will be enforced on block diagonal
    for blk = 1:obj.Nt
        I = (blk-1)*obj.Nr + 1:blk*obj.Nr; % indices for this block
        A(I,:) = obj.enforceBlockIC(blk,A(I,:));
    end
    obj.A = A; obj.B = B;

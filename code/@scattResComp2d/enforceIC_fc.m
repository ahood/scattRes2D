function enforceIC_fc(obj)
    A_fc = obj.A_fc; B_fc = obj.B_fc;
    
    % zero out all rows where interface conditions will be placed
    A_fc = reshape(A_fc,obj.Nr,[]); 
    B_fc = reshape(B_fc,obj.Nr,[]); % put all blocks in one row
    A_fc(obj.valRows,:) = 0; 
    B_fc(obj.valRows,:) = 0; % zero out continuity condition rows
    A_fc(obj.derRows,:) = 0; 
    B_fc(obj.derRows,:) = 0; %          differentiability
    A_fc = reshape(A_fc,obj.Nr*obj.Nt,[]); 
    B_fc = reshape(B_fc,obj.Nr*obj.Nt,[]);
    
    % interface conditions will be enforced on block diagonal
    for blk = 1:obj.Nt
        I = (blk-1)*obj.Nr + 1:blk*obj.Nr; % indices for this block
        A_fc(I,I) = obj.enforceBlockIC_fc(blk,A_fc(I,I));
    end
    obj.A_fc = A_fc; obj.B_fc = B_fc;

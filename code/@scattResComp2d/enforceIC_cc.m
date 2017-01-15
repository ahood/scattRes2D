function enforceIC_cc(obj)
    A_cc = obj.A_cc; B_cc = obj.B_cc;
    
    % zero out all rows where interface conditions will be placed
    A_cc = reshape(A_cc,obj.Nr,[]); 
    B_cc = reshape(B_cc,obj.Nr,[]); % put all blocks in one row
    A_cc(obj.valRows_cc,:) = 0; 
    B_cc(obj.valRows_cc,:) = 0; % zero out continuity condition rows
    A_cc(obj.derRows_cc,:) = 0; 
    B_cc(obj.derRows_cc,:) = 0; %          differentiability
    A_cc = reshape(A_cc,obj.Nr*obj.Nt,[]); 
    B_cc = reshape(B_cc,obj.Nr*obj.Nt,[]);
    
    % interface conditions will be enforced on block diagonal
    for blk = 1:obj.Nt
        I = (blk-1)*obj.Nr + 1:blk*obj.Nr; % indices for this block
        A_cc(I,I) = obj.enforceBlockIC_cc(blk,A_cc(I,I));
    end
    obj.A_cc = A_cc; obj.B_cc = B_cc;

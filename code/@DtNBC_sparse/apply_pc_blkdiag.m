function x = apply_pc_blkdiag(obj,x0,k)
% performs M\x0 for x0 a vector and M the block diagonal of T(k)

if k ~= obj.blk_precond_k % have to remake the blocks before applying
    obj.blk_precond_k = k;
    dtncoeffs = diag( obj.DtNcoeffs(obj.Ns,k,obj.r(end)) );
    C = sparse(obj.U*dtncoeffs*obj.Uinv); 
    c = diag(C);

    for j = 1:obj.Nt
        Mj = obj.A_blks_common + diag(obj.Vhat(:,j) - k^2);

        % interface conditions (if any) -- THIS IS PROBABLY WRONG FOR MORE THAN
        % ONE INTERFACE
        for reg = 2:length(obj.Nrs)
            row = obj.Nrs(reg-1);
            Mj(row,:) = 0;
            Mj(row,row:row+1) = [-1,1];
            Mj(row+1,     1:row) =  obj.Dr_I(row  ,     1:row);
            Mj(row+1, row+1:end) = -obj.Dr_I(row+1, row+1:end);
        end

        % boundary condition
        Mj(end,:) = obj.Dr_I(end,:);
        Mj(end,end) = Mj(end,end) - c(j); % boundary condition
        [L,U] = lu(Mj);
        obj.M_blks{j} = struct('L',L,'U',U);
    end
end

xhat = reshape(x0,obj.Nr,obj.Nt);
for j = 1:obj.Nt
    L = obj.M_blks{j}.L;
    U = obj.M_blks{j}.U;
    xhat(:,j) = U\(L\xhat(:,j));   % do Xj\xhat(:,j)
end
x = xhat(:);


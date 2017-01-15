function y = apply_pc_blkdiag(obj,x0,k)
% performs M\x0 for x0 a vector and M the block diagonal of Tfull(k)

if k ~= obj.blk_precond_k % have to remake the blocks before applying
    obj.blk_precond_k = k;

    for j = 1:obj.dtns.Nt
        M11j = obj.dtns.A_blks_common + diag(obj.dtns.Vhat(:,j) - k^2);

        % interface conditions (if any) -- THIS IS PROBABLY WRONG FOR MORE THAN
        % ONE INTERFACE
        for reg = 2:length(obj.dtns.Nrs)
            row = obj.dtns.Nrs(reg-1);
            M11j(row,:) = 0;
            M11j(row,row:row+1) = [-1,1];
            M11j(row+1,     1:row) =  obj.dtns.Dr_I(row  ,     1:row);
            M11j(row+1, row+1:end) = -obj.dtns.Dr_I(row+1, row+1:end);
        end

        % derivative part of boundary condition
        M11j(end,:) = obj.dtns.Dr_I(end,:);

        % get LU factors for later
        [L,U] = lu(M11j);
        obj.M11_blks{j} = struct('L',L,'U',U);
    end
end

% split into two pieces
x1 = x0(obj.sc.I1);
x2 = x0(obj.sc.I2);

% apply M11 to x1
x1hat = reshape(x1,obj.dtns.Nr,obj.dtns.Nt);
for j = 1:obj.dtns.Nt
    L = obj.M11_blks{j}.L;
    U = obj.M11_blks{j}.U;
    x1hat(:,j) = U\(L\x1hat(:,j));   % do M11j\x1hat(:,j)
end
y1 = x1hat(:);
y2 = x2./(obj.A22diag - k^2);

y = [y1; y2];


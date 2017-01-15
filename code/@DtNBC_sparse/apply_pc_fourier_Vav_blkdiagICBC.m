function y = apply_pc_fourier_Vav_blkdiagICBC(obj,x0,k,kind)
% performs M\x0 for x0 a vector and M comes from approximating VvaluesVec
% as I kron Vav, where Vav is the average of the potential on concentric
% circles, transforming to fourier space, acting in fourier space, and
% transforming back. The action in fourier space is block diagonal due to
% the approximation to the potential.
%
% If kind = 'mult' we multiply by the approx to T(k),
% if kind = 'inv'  we multiply by (approx to T(k))^{-1}.

y = obj.apply_Uinv_kron_I(x0);
yhat = reshape(y, obj.Nr, obj.Nt);

if strcmp(kind,'mult')
    rowChange = obj.ArowChange - k^2*obj.BrowChange;
    eN = 0*yhat(:,1); eN(end) = 1;
    dtncoeffs = obj.DtNcoeffs(obj.Ns,k,obj.r(end));

    for j = 1:obj.Nt
        yj = yhat(:,j);
        Uinvj_LP = kron(obj.Uinv(j,:), obj.localPlacement);
        Uj_I     = kron(obj.U(:,j)   , eye(obj.Nr));

        yhat(:,j) = -(obj.apply_Drr_I(yj) + obj.apply_Dr_I(yj)./obj.r) - ...
                     (obj.apply_Drr_J(yj) + obj.apply_Dr_J(yj)./obj.r)*(-1)^obj.Ns(j) + ...
                    obj.Ns(j)^2*yj./obj.r.^2 + ...
                    (obj.Vav - k^2).*yj + ...
                    Uinvj_LP*( rowChange*(Uj_I*yj) ) - ...
                    dtncoeffs(j)*eN.*yj;
    end
end

if strcmp(kind,'inv')
    % do LU factorizations for all the blocks
    if k ~= obj.blk_fourier_precond_k % have to remake the blocks before applying
        obj.blk_fourier_precond_k = k;
        dtncoeffs = obj.DtNcoeffs(obj.Ns,k,obj.r(end));

        Mshift = obj.ArowChange - k^2*obj.BrowChange;
        % make each block
        for j = 1:obj.Nt
            
            % Do Xj\yhat(:,j) where
            % Xj = -(Drr_I + R*Dr_I)
            %      -(Drr_J + R*Dr_J)*(-1)^(Ns(j))
            %      +Ns(j)^2*R^2 + diag(Vav-k^2) - dtncoeffs(j)*eN*eN^T
            %      +block diag part of row changes in Fourier space

            Xj = -(obj.Drr_I + obj.R*obj.Dr_I) ...
                 -(obj.Drr_J + obj.R*obj.Dr_J)*(-1)^(obj.Ns(j)) ...
                 +obj.Ns(j)^2*obj.R^2 + diag(obj.Vav - k^2);

            % IC and derivative part of scattering boundary condition on block diag
            % DSB: This is still slower than it should be -- the kron is
            %      a slow way to handle the scaling.
            Uinvj_LP = kron(obj.Uinv(j,:), obj.localPlacement);
            Uj_I     = kron(obj.U(:,j)   , speye(obj.Nr));
            Xj = Xj + Uinvj_LP*Mshift*Uj_I;
            
            % DtN map part of scattering boundary condition
            Xj(end,end) = Xj(end,end) - dtncoeffs(j);

            % do LU and save it
            [L,U] = lu(Xj);
            obj.fourierM_blks{j} = struct('L',L,'U',U);
        end
    end

    for j = 1:obj.Nt
        % Do Xj\yhat(:,j) where
        % Xj = -(Drr_I + R*Dr_I)
        %      -(Drr_J + R*Dr_J)*(-1)^(Ns(j))
        %      +Ns(j)^2*R^2 + diag(Vav-k^2) - dtncoeffs(j)*eN*eN^T.

        L = obj.fourierM_blks{j}.L;
        U = obj.fourierM_blks{j}.U;
        yhat(:,j) = U\(L\(yhat(:,j)));
    end
end
y = obj.apply_U_kron_I(yhat(:));

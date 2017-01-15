function y = apply_scpc_fourier_Vav_blkdiagICBC(obj,x0,k,kind)
% performs M\x0 for x0 a vector and M comes from approximating VvaluesVec
% as I kron Vav, where Vav is the average of the potential on concentric
% circles, transforming to fourier space, acting in fourier space, and
% transforming back. The action in fourier space is block diagonal due to
% the approximation to the potential.
%
% If kind = 'mult' we multiply by the approx to T(k), 
% if kind = 'inv'  we multiply by (approx to T(k))^{-1}.

y = obj.dtns.apply_Uinv_kron_I(x0);
yhat = reshape(y, obj.dtns.Nr, obj.dtns.Nt);

if strcmp(kind,'mult')
    rowChange = obj.dtns.ArowChange - k^2*obj.dtns.BrowChange;
    eN = 0*yhat(:,1); eN(end) = 1;
    
    for j = 1:obj.dtns.Nt
        yj = yhat(:,j);
        Uinvj_LP = kron(obj.dtns.Uinv(j,:), obj.dtns.localPlacement);
        Uj_I     = kron(obj.dtns.U(:,j)   , eye(obj.dtns.Nr));
        approx_dtncoeffj = obj.ratf(j).eval_at(k^2);
        
        yhat(:,j) = -(obj.dtns.apply_Drr_I(yj) + obj.dtns.apply_Dr_I(yj)./obj.dtns.r) - ...
                     (obj.dtns.apply_Drr_J(yj) + obj.dtns.apply_Dr_J(yj)./obj.dtns.r)*(-1)^obj.dtns.Ns(j) + ...
                    obj.dtns.Ns(j)^2*yj./obj.dtns.r.^2 + ...
                    (obj.dtns.Vav - k^2).*yj + ...
                    Uinvj_LP*( rowChange*(Uj_I*yj) ) - ...
                    approx_dtncoeffj*eN.*yj;
    end
end

if strcmp(kind,'inv')
    % do LU factorizations for all the blocks
    if k ~= obj.sc_blk_fourier_precond_k % have to remake the blocks before applying
        obj.sc_blk_fourier_precond_k = k;

        % make each block
        for j = 1:obj.dtns.Nt
            % Do Xj\yhat(:,j) where
            % Xj = -(Drr_I + R*Dr_I)
            %      -(Drr_J + R*Dr_J)*(-1)^(Ns(j))
            %      +Ns(j)^2*R^2 + diag(Vav-k^2) - dtncoeffs(j)*eN*eN^T 
            %      +block diag part of row changes in Fourier space

            Xj = -(obj.dtns.Drr_I + obj.dtns.R*obj.dtns.Dr_I) ...
                 -(obj.dtns.Drr_J + obj.dtns.R*obj.dtns.Dr_J)*(-1)^(obj.dtns.Ns(j)) ...
                 +obj.dtns.Ns(j)^2*obj.dtns.R^2 + diag(obj.dtns.Vav - k^2);

            % IC and derivative part of scattering boundary condition on block diag
            Uinvj_LP = kron(obj.dtns.Uinv(j,:), obj.dtns.localPlacement);
            Uj_I     = kron(obj.dtns.U(:,j)   , eye(obj.dtns.Nr));
            Xj = Xj + Uinvj_LP*(obj.dtns.ArowChange - k^2*obj.dtns.BrowChange)*Uj_I;

            % DtN map part of scattering boundary condition
            approx_dtncoeffj = obj.ratf(j).eval_at(k^2);
            Xj(end,end) = Xj(end,end) - approx_dtncoeffj;

            % do LU and save it
            [L,U] = lu(Xj);
            obj.sc_fourierM_blks{j} = struct('L',L,'U',U);
        end
    end

    for j = 1:obj.dtns.Nt
        % Do Xj\yhat(:,j) where
        % Xj = -(Drr_I + R*Dr_I)
        %      -(Drr_J + R*Dr_J)*(-1)^(Ns(j))
        %      +Ns(j)^2*R^2 + diag(Vav-k^2) - dtncoeffs(j)*eN*eN^T.

        L = obj.sc_fourierM_blks{j}.L;
        U = obj.sc_fourierM_blks{j}.U;
        yhat(:,j) = U\(L\yhat(:,j));
    end 
end
y = obj.dtns.apply_U_kron_I(yhat(:));

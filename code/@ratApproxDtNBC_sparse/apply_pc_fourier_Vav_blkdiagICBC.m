function w = apply_pc_fourier_Vav_blkdiagICBC(obj,x0,k,kind)
% performs M\x0 for x0 a vector and 
% M = [A11 - k^2*B11, A12; A21, diag(A22diag-k^2)] using
% A11 - k^2*B11 = kron(U,I)*blkdiag(Y1,...,YN)kron(Uinv,I) 
% (same as blkdiag(X1,...,XN) from
% version in DtNBC_sparse but without DtN map term).
% Inverting M with 2x2 block formula
% inv[A,B;C,D] = [inv(S), -S\(B*inv(D)); etc.].
%
% If kind = 'mult' we multiply by M
% if kind = 'inv'  we multiply by M^{-1}.

x1 = x0(obj.sc.I1);
x2 = x0(obj.sc.I2);

% apply [ kron(Uinv,I), 0; 0, I]
y1 = obj.dtns.apply_Uinv_kron_I(x1);
y2 = x2;

% need this for applying block diagonal stuff
y1hat = reshape(y1, obj.dtns.Nr, obj.dtns.Nt);

if strcmp(kind,'mult')
    rowChange = obj.dtns.ArowChange - k^2*obj.dtns.BrowChange;

    % apply Y = blkdiag(Y1,...,YN) to y1
    Yy1 = y1hat;
    for j = 1:obj.dtns.Nt
        y1j = y1hat(:,j);
        Uinvj_LP = kron(obj.dtns.Uinv(j,:), obj.dtns.localPlacement);
        Uj_I     = kron(obj.dtns.U(:,j)   , eye(obj.dtns.Nr));
        
        Yy1(:,j) = -(obj.dtns.apply_Drr_I(y1j) + obj.dtns.apply_Dr_I(y1j)./obj.dtns.r) - ...
                    (obj.dtns.apply_Drr_J(y1j) + obj.dtns.apply_Dr_J(y1j)./obj.dtns.r)*(-1)^obj.dtns.Ns(j) + ...
                    obj.dtns.Ns(j)^2*y1j./obj.dtns.r.^2 + ...
                    (obj.dtns.Vav - k^2).*y1j + ...
                    Uinvj_LP*( rowChange*(Uj_I*y1j) );
    end
    z1 = Yy1(:) + obj.apply_fourierA12(y2);
    z2 = obj.apply_fourierA21(y1) + (obj.A22diag - k^2).*y2;
end

if strcmp(kind,'inv')
    % do LU factorizations for all the blocks of fourierS, where
    % S = kron(U,I)*Y*kron(Uinv,I) - A12*inv(D22-k^2*I)*A21 and
    % fourierS = Y - A12orig*inv(D22-k^2*I)*A21orig
    if k ~= obj.blk_fourier_precond_k % have to remake the blocks before applying
        obj.blk_fourier_precond_k = k;

        for j = 1:obj.dtns.Nt
            % fourierSj = Yj - block j of (fourierA12*inv(D22-k^2*I)*fourierA21)
            % where
            % Yj = -(Drr_I + R*Dr_I)
            %      -(Drr_J + R*Dr_J)*(-1)^(Ns(j))
            %      +Ns(j)^2*R^2 + diag(Vav-k^2)
            %      +block diag part of row changes in Fourier space

            % make blocks of Y
            Yj = -(obj.dtns.Drr_I + obj.dtns.R*obj.dtns.Dr_I) ...
                 -(obj.dtns.Drr_J + obj.dtns.R*obj.dtns.Dr_J)*(-1)^(obj.dtns.Ns(j)) ...
                 +obj.dtns.Ns(j)^2*obj.dtns.R^2 + diag(obj.dtns.Vav - k^2);

            % IC and derivative part of scattering boundary condition on block diag
            Uinvj_LP = kron(obj.dtns.Uinv(j,:), obj.dtns.localPlacement);
            Uj_I     = kron(obj.dtns.U(:,j)   , eye(obj.dtns.Nr));
            Yj = Yj + Uinvj_LP*(obj.dtns.ArowChange - k^2*obj.dtns.BrowChange)*Uj_I;

            % DtN map part of scattering boundary condition
            fourierSj = Yj;
%             approx_dtncoeffj = sum(obj.ratf(j).w./(obj.ratf(j).z - k^2));
            approx_dtncoeffj = obj.ratf(j).eval_at(k^2);
            fourierSj(end,end) = fourierSj(end,end) - approx_dtncoeffj;
            
            % do LU and save it
            [L,U] = lu(fourierSj);
            obj.fourierS_blks{j} = struct('L',L,'U',U);
        end
    end

    invSy1hat = y1hat;
    A12inv22y2 = obj.apply_fourierA12(y2./(obj.A22diag - k^2));
    A12inv22y2hat = reshape(A12inv22y2, obj.dtns.Nr, obj.dtns.Nt);
    invSA12inv22y2hat = A12inv22y2hat;
    for j = 1:obj.dtns.Nt
        % Do [ inv(fourierS)                          ,   -inv(fourierS)*fourierA12*inv(D22 - k^2*I); 
        %     -inv(D22-k^2*I)*fourierA21*inv(fourierS), inv(D22-k^2*I) +inv(D22-k^2*I)*fourierA21*inv(fourierS)*fourierA12*inv(D22-k^2*I)]
        % times [y1; y2] given LU factorization of fourierS.

        L = obj.fourierS_blks{j}.L;
        U = obj.fourierS_blks{j}.U;
        
        invSy1hat(:,j) = U\(L\y1hat(:,j)); % fourierS\y1
        invSA12inv22y2hat(:,j) = U\(L\A12inv22y2hat(:,j)); % fourierS\(fourierA12*inv(D22-k^2*I)*y2)
    end 
    invSy1 = invSy1hat(:);
    invSA12inv22y2 = invSA12inv22y2hat(:);
    z1 = invSy1 - invSA12inv22y2;
    z2 = (-obj.apply_fourierA21(invSy1) + y2 + obj.apply_fourierA21(invSA12inv22y2))./(obj.A22diag - k^2);
end
% apply [kron(U,I), 0; 0, I]
w = [obj.dtns.apply_U_kron_I(z1); z2];

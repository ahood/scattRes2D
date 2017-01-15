function [x,DrIxhat,DrJxhat] = apply_op_no_bcs(obj,x0,E0,VvaluesVec)
% Do x = (A - E0*B)*x0 where no boundary conditions have been applied.
% B is just the identity.

    % need this as input to apply_kron
    xhat = reshape(x0, obj.Nr, []); % for applying kronecker products

    x = VvaluesVec.*x0; % potential term
    x = x - obj.apply_kron(obj.Dtt,diag(1./obj.r.^2),xhat); % Dtt term
    x = x - obj.apply_J_kron_Drr_J(xhat); % Drr_J term
    x = x - obj.apply_I_kron_Drr_I(xhat); % Drr_I term
    DrIxhat = obj.apply_I_kron_Dr_I(xhat);
    DrJxhat = obj.apply_J_kron_Dr_J(xhat);
    x = x - (DrIxhat + DrJxhat)./obj.rs;

    x = x - E0*x0;
    
    % apply C0 condition at interfaces
    x = reshape(x, obj.Nr, []);
%     for valRow = obj.valRows
    i1s = min(obj.valRows,obj.derRows);
    i2s = max(obj.valRows,obj.derRows);
    for j = 1:length(obj.valRows)
        valRow = obj.valRows(j);
        i1 = i1s(j); i2 = i2s(j);
%         x(valRow,:) = xhat(valRow+1,:) - xhat(valRow,:);
        x(valRow,:) = xhat(i2,:) - xhat(i1,:);
    end
    x = x(:);

    % apply C1 condition at interfaces
%     for derRow = obj.derRows
    for j = 1:length(obj.derRows)
        derRow = obj.derRows(j);
        for ii = 1:obj.Nt
            row = derRow + (ii-1)*obj.Nr; % this row of x is getting replaced
            i1 = i1s(j) + (ii-1)*obj.Nr;
            i2 = i2s(j) + (ii-1)*obj.Nr;
%             x(row) = DrIxhat(row-1) - DrIxhat(row) + ...
%                      DrJxhat(row-1) - DrJxhat(row);
            x(row) = DrIxhat(i1) - DrIxhat(i2) + ...
                     DrJxhat(i1) - DrJxhat(i2);
        end
    end

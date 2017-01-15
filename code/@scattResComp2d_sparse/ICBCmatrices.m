function ICBCmatrices(obj)
% Computes ArowChange, BrowChange such that
% A = Aorig + kron(I,obj.localPlacement)*ArowChange and
% B = Borig + kron(I,obj.localPlacement)*BrowChange

% compute ArowChange
It = speye(obj.Nt);
Ithalf = speye(obj.Nt/2);
Jt = [0*Ithalf,   Ithalf; ...
        Ithalf, 0*Ithalf];

n = length(obj.VvaluesVec);
V = spdiags(obj.VvaluesVec,0,n,n);
    
Aorig_rows = -(kron(It,obj.R(obj.localICBCrows,:)*obj.Dr_I + obj.Drr_I(obj.localICBCrows,:)) + ...
               kron(Jt,obj.R(obj.localICBCrows,:)*obj.Dr_J + obj.Drr_J(obj.localICBCrows,:)) + ...
               kron(obj.Dtt,obj.R(obj.localICBCrows,:)*obj.R)) + ...
             V(obj.globalICBCrows,:);
         
A_rows = zeros(obj.Nt*(2*length(obj.valRows) + 1), obj.Nr*obj.Nt);
Dr_rows = kron(It,obj.Dr_I(obj.localICBCrows,:)) + ...
          kron(Jt,obj.Dr_J(obj.localICBCrows,:));

for blk = 1:obj.Nt
    row_shift = (blk-1)*(2*length(obj.valRows) + 1);
    col_shift = (blk-1)*obj.Nr;
    % interface conditions  
    for j = 1:length(obj.valRows) % interface j
%         C0_row = row_shift + 2*j-1; C1_row = row_shift + 2*j;
%         vj = col_shift + obj.valRows(j); dj = col_shift + obj.derRows(j);
%         A_rows(C0_row,[vj,dj]) = [-1,1]; % continuity
%         A_rows(C1_row,:) = Dr_rows(C0_row,:) - Dr_rows(C1_row,:); % C1
        vj = obj.valRows(j); dj = obj.derRows(j);
        if vj < dj
            C0_row = row_shift + 2*j-1; C1_row = row_shift + 2*j;
        else
            C1_row = row_shift + 2*j-1; C0_row = row_shift + 2*j;
        end
        i1 = min(col_shift + vj, col_shift + dj);
        i2 = max(col_shift + vj, col_shift + dj);
        A_rows(C0_row,[i1,i2]) = [-1,1]; % continuity
        A_rows(C1_row,:) = Dr_rows(min(C0_row,C1_row),:) - ...
                           Dr_rows(max(C0_row,C1_row),:); % C1
    end
    % boundary condition
    A_rows(row_shift + 2*j+1,:) = Dr_rows(row_shift + 2*j+1,:);
end

obj.ArowChange = A_rows - Aorig_rows;

% compute BrowChange
Ir = speye(obj.Nr);
Borig_rows = kron(It,Ir(obj.localICBCrows,:));

B_rows = 0*A_rows;

obj.BrowChange = B_rows - Borig_rows;

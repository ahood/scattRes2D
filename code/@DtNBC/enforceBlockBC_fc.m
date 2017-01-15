function X = enforceBlockBC_fc(obj,blk,X)
% Takes block row for certain angle and enforces derivative part of
% boundary condition.

    % select a differentiation matrix
    if mod(obj.Nt/2,2) == 0 % last one is even
        if mod(blk,2) == 0, D = obj.Drs_even;
        else                D = obj.Drs_odd;
        end
    else
        if mod(blk,2) == 0, D = obj.Drs_odd;
        else                D = obj.Drs_even;
        end
    end

    % derivative part of boundary condition
    X(end,:) = D(end,:);

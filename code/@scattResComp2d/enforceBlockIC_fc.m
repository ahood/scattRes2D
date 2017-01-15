function X = enforceBlockIC_fc(obj,blk,X)
% Takes block row for certain angle and enforces interface conditions.

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

    % interface conditions    
    for j = 1:length(obj.valRows) % interface j
        vj = obj.valRows(j); dj = obj.derRows(j);
        i1 = min(vj,dj); i2 = max(vj,dj);
        X(vj,[i1,i2]) = [-1,1]; % continuity
        X(dj,:) = D(i1,:) - D(i2,:); % C1

%         X(vj,[vj,dj]) = [-1,1]; % continuity
%         X(dj,:) = D(vj,:) - D(dj,:);
    end

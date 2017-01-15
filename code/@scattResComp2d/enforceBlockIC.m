function X = enforceBlockIC(obj,blk,X)
% Takes block row for certain angle and enforces interface conditions.

    shift = (blk-1)*obj.Nr;

    % interface conditions    
    for j = 1:length(obj.valRows) % interface j
        vj = obj.valRows(j); dj = obj.derRows(j);
%         X(vj,shift + [vj,dj]) = [-1,1]; % continuity
%         X(dj,:) = obj.Dr(shift+vj,:) - obj.Dr(shift+dj,:); % C1
        i1 = min(vj,dj); i2 = max(vj,dj);
        X(vj,shift + [i1,i2]) = [-1,1]; % continuity
        X(dj,:) = obj.Dr(shift+i1,:) - obj.Dr(shift+i2,:); % C1
    end

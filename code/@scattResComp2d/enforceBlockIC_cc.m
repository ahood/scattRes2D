function X = enforceBlockIC_cc(obj,blk,X)
% Takes block row for certain angle and enforces interface conditions.

    shift = (blk-1)*obj.Nr;
    
    % check if this block corresponds to an even or odd function
    if mod(obj.Nt/2,2) == 0 % last block is even
        if mod(blk,2) == 0 
            Winv_blk = obj.Winv_even; 
            Drs_blk  = obj.Drs_even;
        else
            Winv_blk = obj.Winv_odd;
            Drs_blk  = obj.Drs_odd;
        end
    else
        if mod(blk,2) == 0
            Winv_blk = obj.Winv_odd;
            Drs_blk  = obj.Drs_odd;
        else
            Winv_blk = obj.Winv_even;
            Drs_blk = obj.Drs_even;
        end
    end

    % interface conditions (imposed in Fourier space)
    for j = 1:length(obj.valRows_cc) % interface j
        % on jth interface
        rowL = obj.valRows_cc(j);     % last point in interval to left
        rowR = obj.valRows_cc(j) + 1; %                           right
        X(obj.valRows_cc(j),:) = ...
            Winv_blk(rowR,:) - ...
            Winv_blk(rowL,:);

        X(obj.derRows_cc(j),:) = ...
            Drs_blk(rowR,:)*Winv_blk - ...
            Drs_blk(rowL,:)*Winv_blk;
    end

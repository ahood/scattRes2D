function X = enforceBlockBC_cc(obj,blk,X)
% Takes block row for certain angle and enforces derivative part of
% boundary condition.

    % select a differentiation matrix mapping cheb coeffs to derivative in 
    % fourier basis
    if mod(obj.Nt/2,2) == 0 % last one is even
        if mod(blk,2) == 0
            X(end,:) = obj.Drs_even(end,:)*obj.Winv_even;
        else                
            X(end,:) = obj.Drs_odd(end,:)*obj.Winv_odd;
        end
    else
        if mod(blk,2) == 0
            X(end,:) = obj.Drs_odd(end,:)*obj.Winv_odd;
        else                
            X(end,:) = obj.Drs_even(end,:)*obj.Winv_even;
        end
    end

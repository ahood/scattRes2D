function Y = apply_Winv_on_right(obj,X)
% Let W be the map from values to cheb coeffs (applied using Fourier
% transform and symmetry considerations). 
% Then Winv = blkdiag(Winv_even, Winv_odd , ...) or
%      Winv = blkdiag(Winv_odd , Winv_even, ...),
% depending on the parity of obj.Nt/2.
% Returns X*Winv.

% ---------------------------------
% TODO
% Shouldn't deriv also be an argument?
% ---------------------------------

    XT = X.';
    XT_reshaped = reshape(XT, 2*obj.Nr, obj.Nr*obj.Nt*obj.Nt/2);
    X_reshaped = XT_reshaped.';
    if mod(obj.Nt/2,2) == 0
        Y_reshaped = X_reshaped*blkdiag(obj.Winv_odd , obj.Winv_even);
    else
        Y_reshaped = X_reshaped*blkdiag(obj.Winv_even, obj.Winv_odd );
    end
    YT_reshaped = Y_reshaped.';
    YT = reshape(YT_reshaped, obj.Nr*obj.Nt, obj.Nr*obj.Nt);
    Y = YT.';

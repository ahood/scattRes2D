function enforceBC_full_cc(obj)
% Enforce Dirichlet boundary conditions in Chebyshev basis.

    % put last row of Winv blocks in A_cc
    if mod(obj.Nt/2,2) == 0 % last block is even
        Winv1 = obj.Winv_odd;  Winv2 = obj.Winv_even;
    else
        Winv1 = obj.Winv_even; Winv2 = obj.Winv_odd;
    end

    for j = 1:2:obj.Nt
        obj.Afull_cc_blks{j  }(end,:) = Winv1(end,:);
        obj.Afull_cc_blks{j+1}(end,:) = Winv2(end,:);
    end
    for j = 1:obj.Nt
        obj.Bfull_cc_blks{j}(end,:) = 0;
    end
        

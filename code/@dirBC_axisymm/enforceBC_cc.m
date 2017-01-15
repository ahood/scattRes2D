function enforceBC_cc(obj)
% Enforce Dirichlet boundary conditions in Chebyshev basis.
    
    if mod(obj.Nt/2,2) == 0 % last block is even
        Winv1 = obj.Winv_odd;  Winv2 = obj.Winv_even;
    else
        Winv1 = obj.Winv_even; Winv2 = obj.Winv_odd;
    end

    % put last row of Winv blocks in A_cc
    for j = 1:2:obj.Nt-1
        obj.A_cc_blks{j  }(end,:) = Winv1(end,:);
        obj.A_cc_blks{j+1}(end,:) = Winv2(end,:);
    end
    for j = 1:obj.Nt
        obj.B_cc_blks{j  }(end,:) = 0;
    end

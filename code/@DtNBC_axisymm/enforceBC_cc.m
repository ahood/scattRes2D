function enforceBC_cc(obj)
% enforces boundary conditions in cheb basis

    if mod(obj.Nt/2,2) == 0 % last block even
        D1 = obj.Drs_odd;  Winv1 = obj.Winv_odd;  
        D2 = obj.Drs_even; Winv2 = obj.Winv_even;
    else
        D1 = obj.Drs_even; Winv1 = obj.Winv_even; 
        D2 = obj.Drs_odd;  Winv2 = obj.Winv_odd;
    end

    for j = 1:2:obj.Nt-1
%         obj.A_cc_blks{j  }(end,:) = D1(end,:)*Winv1;
%         obj.A_cc_blks{j+1}(end,:) = D2(end,:)*Winv2;
        obj.A_cc_blks{j  }(end,:) = -D1(end,:)*Winv1;
        obj.A_cc_blks{j+1}(end,:) = -D2(end,:)*Winv2;
    end
    for j = 1:obj.Nt
        obj.B_cc_blks{j}(end,:) = 0;
    end

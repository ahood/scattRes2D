function C_cc(obj,k)

    d = obj.DtNcoeffs(obj.Ns,k,obj.r(end));
    
    if mod(obj.Nt/2,2) == 0 % last block even
        Winv1 = obj.Winv_odd;  Winv2 = obj.Winv_even;
    else
        Winv1 = obj.Winv_even; Winv2 = obj.Winv_odd;
    end
    
    for j = 1:2:obj.Nt-1
        obj.C_cc_blks{j  } = spalloc(obj.Nr,obj.Nr,obj.Nr);
%         obj.C_cc_blks{j  }(end,:) = -d(j  )*Winv1(end,:);
        obj.C_cc_blks{j  }(end,:) = d(j  )*Winv1(end,:);

        obj.C_cc_blks{j+1} = spalloc(obj.Nr,obj.Nr,obj.Nr);
%         obj.C_cc_blks{j+1}(end,:) = -d(j+1)*Winv2(end,:);
        obj.C_cc_blks{j+1}(end,:) = d(j+1)*Winv2(end,:);
    end
end

function enforceIC_full_cc(obj)

    if mod(obj.Nt/2,2) == 0 % last block is even
        Winv1 = obj.Winv_odd;  Winv2 = obj.Winv_even;
        D1 = obj.Drs_odd;      D2 = obj.Drs_even;
    else
        Winv1 = obj.Winv_even; Winv2 = obj.Winv_odd;
        D1 = obj.Drs_even;     D2 = obj.Drs_odd;
    end        

    % make the continuity and derivative rows to plug in
    cont1 = zeros(length(obj.valRows_cc), obj.Nr);
    cont2 = cont1;
    deriv1 = cont1;
    deriv2 = cont1;
    for j = 1:length(obj.valRows_cc)
        vj = obj.valRows_cc(j); dj = obj.derRows_cc(j);
        cont1(j,:) = Winv1(vj+1,:) - Winv1(vj,:);
        cont2(j,:) = Winv2(vj+1,:) - Winv2(vj,:);
        deriv1(j,:) = (D1(vj+1,:) - D1(vj,:))*Winv1;
        deriv2(j,:) = (D2(vj+1,:) - D2(vj,:))*Winv2;
    end
    
    for j = 1:2:obj.Nt-1
        obj.Afull_cc_blks{j  }(obj.valRows_cc,:) = cont1;
        obj.Afull_cc_blks{j  }(obj.derRows_cc,:) = deriv1;
        obj.Afull_cc_blks{j+1}(obj.valRows_cc,:) = cont2;
        obj.Afull_cc_blks{j+1}(obj.derRows_cc,:) = deriv2;
    end
    for j = 1:obj.Nt
        obj.Bfull_cc_blks{j}([obj.valRows_cc(:); obj.derRows_cc(:)],:) = 0;
    end

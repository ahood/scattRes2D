function enforceIC_full_fc(obj)

    if mod(obj.Nt/2,2) == 0 % last one is even
        D1 = obj.Drs_odd;  D2 = obj.Drs_even;
    else
        D1 = obj.Drs_even; D2 = obj.Drs_odd;
    end

    % make the continuity and derivative rows to plug in
    continuity = zeros(length(obj.valRows), obj.Nr);
    deriv1 = continuity;
    deriv2 = continuity;
    for j = 1:length(obj.valRows)
        vj = obj.valRows(j); dj = obj.derRows(j);
        continuity(j,[vj,dj]) = [-1,1];
        deriv1(j,:) = D1(vj,:) - D1(dj,:);
        deriv2(j,:) = D2(vj,:) - D2(dj,:);
    end
    
    for j = 1:2:obj.Nt-1
        obj.Afull_fc_blks{j  }(obj.valRows,:) = continuity;
        obj.Afull_fc_blks{j  }(obj.derRows,:) = deriv1;
        obj.Afull_fc_blks{j+1}(obj.valRows,:) = continuity;
        obj.Afull_fc_blks{j+1}(obj.derRows,:) = deriv2;
    end
    for j = 1:obj.Nt
        obj.Bfull_fc_blks{j}([obj.valRows(:); obj.derRows(:)],:) = 0;
    end

function enforceBC_fc(obj)

    if mod(obj.Nt/2,2) == 0 % last block even
        D1 = obj.Drs_odd;  D2 = obj.Drs_even;
    else
        D1 = obj.Drs_even; D2 = obj.Drs_odd;
    end

    for j = 1:2:obj.Nt-1
%         obj.A_fc_blks{j  }(end,:) = D1(end,:);
%         obj.A_fc_blks{j+1}(end,:) = D2(end,:);
        obj.A_fc_blks{j  }(end,:) = -D1(end,:);
        obj.A_fc_blks{j+1}(end,:) = -D2(end,:);
    end
    for j = 1:obj.Nt
        obj.B_fc_blks{j}(end,:) = 0;
    end

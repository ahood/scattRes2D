function enforceBC_full_fc(obj)

    for j = 1:obj.Nt
        obj.Afull_fc_blks{j}(end,:) = 0;
        obj.Afull_fc_blks{j}(end,end) = 1;
        obj.Bfull_fc_blks{j}(end,:) = 0;
    end

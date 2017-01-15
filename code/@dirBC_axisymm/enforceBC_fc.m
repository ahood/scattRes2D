function enforceBC_fc(obj)

    for j = 1:obj.Nt
        obj.A_fc_blks{j}(end,:) = 0;
        obj.A_fc_blks{j}(end,end) = 1;
        obj.B_fc_blks{j}(end,:) = 0;
    end

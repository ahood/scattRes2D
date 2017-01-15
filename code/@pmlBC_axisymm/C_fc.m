function C_fc(obj,k)

    for j = 1:obj.Nt
        obj.C_fc_blks{j} = ...
            obj.pieces.update(obj.Afull_fc_blks{j} - k^2*...
                              obj.Bfull_fc_blks{j});
    end
function C_cc(obj,k)

    for j = 1:obj.Nt
        obj.C_cc_blks{j} = ...
            obj.pieces.update(obj.Afull_cc_blks{j} - k^2*...
                              obj.Bfull_cc_blks{j});
    end
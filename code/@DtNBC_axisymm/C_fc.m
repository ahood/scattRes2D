function C_fc(obj,k)

    d = obj.DtNcoeffs(obj.Ns,k,obj.r(end));
    
    for j = 1:obj.Nt
        obj.C_fc_blks{j} = spalloc(obj.Nr,obj.Nr,1);
        obj.C_fc_blks{j}(end,end) = d(j);
    end
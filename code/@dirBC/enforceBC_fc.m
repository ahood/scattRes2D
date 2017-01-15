function enforceBC_fc(obj)
    
    obj.A_fc(obj.BCrows,:) = 0;
    obj.B_fc(obj.BCrows,:) = 0;
    
    for row = obj.BCrows
        obj.A_fc(row,row) = 1;
    end
    
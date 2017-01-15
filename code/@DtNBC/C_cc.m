function Ck = C_cc(obj,k)
% returns the nonlinear part of T_cc
    d = obj.DtNcoeffs(obj.Ns,k,obj.r(end));
    n = obj.Nr*obj.Nt;
    Ck = sparse(n,n);
    
    if mod(obj.Nt/2,2) == 0 % last one is even
        Winv1 = obj.Winv_odd;  Winv2 = obj.Winv_even;
    else
        Winv1 = obj.Winv_even; Winv2 = obj.Winv_odd;
    end

    for j = 1:2:obj.Nt-1
        row  = obj.BCrows(j);
        cols = row - obj.Nr + 1 : row;
        Ck(row,cols) = -d(j)*Winv1(end,:);
        
        row  = obj.BCrows(j+1);
        cols = row - obj.Nr + 1 : row;
        Ck(row,cols) = -d(j+1)*Winv2(end,:);
    end

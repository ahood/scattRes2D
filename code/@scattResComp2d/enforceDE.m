function enforceDE(obj,Lr,Vs,coords,r,s)
% Discretize -Laplacian + V - k^2 as A - k^2 B.

    A = sparse(obj.Nr*obj.Nt,obj.Nr*obj.Nt);
    B = speye(obj.Nr*obj.Nt);

    % set up evaluation of V on mesh, one column per angle
    if isnumeric(Vs)
        VvaluesVec = Vs;
    else
        VvaluesVec = obj.valuesVecFromFunCellArray(Vs,coords);
    end
    
    R = sparse(diag(1./s));
    L = Lr + kron(obj.Dtt,R^2);

    A = -sparse(L);
    I = 1:length(A)+1:numel(A);
    A(I) = A(I) + VvaluesVec.';
    
    obj.A = A;
    obj.B = B;
    obj.Aorig = A; % No boundary conditions yet
    obj.Borig = B;
    
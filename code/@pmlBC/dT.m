function dTk = dT(obj,k)
% derivative of T w.r.t. k
    % T = (A11 - zB11) - A12*( (A22-zB22)\A21 ), z = k^2
    % d/dz inv(A22-zB22) = -(A22-zB22)*(-B22)*(A22-zB22),
    % dT/dz = -B11 - A12*inv(A22-zB22)*B22*inv(A22-zB22)
    Tk = obj.T(k);
    I1 = obj.pieces.I1; I2 = obj.pieces.I2;
    A11 = obj.A(I1,I1); A12 = obj.A(I1,I2);
    A21 = obj.A(I2,I1); A22 = obj.A(I2,I2);
    B11 = obj.B(I1,I1); B22 = obj.B(I2,I2);
    z = k^2;
    dTdz = -B11 - (A12/(A22-z*B22))*(B22/(A22-z*B22))*A21;
    dTk = dTdz*2*k;
end

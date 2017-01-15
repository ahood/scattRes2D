function x = apply_T(obj,x0,k)
% Do x = (A - k^2*B)*x0 and enforce Dirichlet BCs at outer edge.
% Inputs:
%  - obj: A PML instance
%  - x0 : vector we're applying obj.T(k) to quickly
%  - k  : energy parameter, i.e. T(k) = A - k^2*B

% Think of Tfull(k) as [M11, M12; M21; M22], so T(k) being the Schur
% complement is M11 - M12*(M22\M21).

% do kronx = kron(Uinv,I)*x0
x0hat = reshape(x0, obj.dtns.Nr, obj.dtns.Nt);
kronx = obj.dtns.apply_Uinv_kron_I(x0hat);

% do x21 = M21*x0
x21 = 0*obj.A22diag;
for jj = 1:length(obj.dtns.BCrows)
    BCcol = obj.dtns.BCrows(jj);
    rows = obj.ratf_ends(jj) - obj.ratf_lens(jj) + 1 : obj.ratf_ends(jj);
    x21(rows) = kronx(BCcol);
end

% do x22 = M22\x21
x22 = x21./(obj.A22diag - k^2);

% do x12 = M12*x22
x12 = 0*x0;
for jj = 1:length(obj.dtns.BCrows)
    BCrow = obj.dtns.BCrows(jj);
    cols = obj.ratf_ends(jj) - obj.ratf_lens(jj) + 1 : obj.ratf_ends(jj);
    x12(BCrow) = obj.ratf(jj).w*x22(cols);
end

% do x12 = kron(U,I)*x12
x12hat = reshape(x12, obj.dtns.Nr, obj.dtns.Nt);
x12 = obj.dtns.apply_U_kron_I(x12hat);

% do x = x11 - x12 where x11 = M11*x0
[x11,DrIxhat,DrJxhat] = obj.dtns.apply_op_no_bcs(x0,k^2,obj.dtns.VvaluesVec);
x11(obj.dtns.BCrows) = DrIxhat(obj.dtns.BCrows) + DrJxhat(obj.dtns.BCrows);

x = x11 - x12;

function y = apply_Tfull(obj,x0,k)
% Do y = (A - k^2*B)*x0 and enforce Dirichlet BCs at outer edge.
% Inputs:
%  - obj: A PML instance
%  - x0 : vector we're applying obj.T(k) to quickly
%  - k  : energy parameter, i.e. T(k) = A - k^2*B

% split x0
x1 = x0(obj.sc.I1); % potential support
x2 = x0(obj.sc.I2); % PML region

% Split y = [y1; y2] same way.

% Get y1 = (A11 - k^2*B11)*x1 + A12*x2 = y11 + y12
[y11,DrIxhat,DrJxhat] = obj.dtns.apply_op_no_bcs(x1,k^2,obj.dtns.VvaluesVec); % no BC's applied yet
y11(obj.dtns.BCrows) = DrIxhat(obj.dtns.BCrows) + DrJxhat(obj.dtns.BCrows);

y12 = 0*x1; % to get right dimensions (Nr*Nt)
for jj = 1:length(obj.dtns.BCrows)
    BCrow = obj.dtns.BCrows(jj);
    cols = obj.ratf_ends(jj) - obj.ratf_lens(jj) + 1 : obj.ratf_ends(jj);
    y12(BCrow) = obj.ratf(jj).w*x2(cols);
end
y12hat = reshape(y12, obj.dtns.Nr, obj.dtns.Nt);
y12 = obj.dtns.apply_U_kron_I(y12hat);

y1 = y11 + y12;

% Get y2 = A21*x1 + (A22 - k^2*B22)*x2 = y21 + y22
x1hat = reshape(x1, obj.dtns.Nr, obj.dtns.Nt);
kronx1hat = obj.dtns.apply_Uinv_kron_I(x1hat);

y21 = 0*x2; % right dimensions
for jj = 1:length(obj.dtns.BCrows)
    BCcol = obj.dtns.BCrows(jj);
    rows = obj.ratf_ends(jj) - obj.ratf_lens(jj) + 1 : obj.ratf_ends(jj);
%     y2(rows) = y2(rows) + x1(BCcol);
    y21(rows) = kronx1hat(BCcol);
end

y22 = (obj.A22diag - k^2).*x2;

y2 = y21 + y22;

y = [y1;y2];

function x = apply_T(obj,x0,k)
% Do x = (A - k^2*B + C(k))*x0 where C(k) is from DtN map.
% Inputs:
%  - obj: A DtN instance
%  - x0 : vector we're applying obj.T(k) to quickly
%  - k  : frequency parameter, i.e. T(k) = A - k^2*B + C(k)

[x,DrIxhat,DrJxhat] = obj.apply_op_no_bcs(x0,k^2,obj.VvaluesVec); % no BC's applied yet

x(obj.BCrows) = DrIxhat(obj.BCrows) + DrJxhat(obj.BCrows) - obj.apply_DtNmap(x0(obj.BCrows),k);



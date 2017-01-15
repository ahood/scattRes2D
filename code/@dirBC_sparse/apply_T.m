function x = apply_T(obj,x0,k)
% Do x = (A - k^2*B)*x0.
% Inputs:
%  - obj: A DtN instance
%  - x0 : vector we're applying obj.T(k) to quickly
%  - k  : frequency parameter, i.e. T(k) = A - k^2*B + C(k)

x = obj.apply_op_no_bcs(x0,k^2,obj.VvaluesVec); % no BC's applied yet

% Dirichlet boundary condition row actions
x(obj.BCrows) = x0(obj.BCrows);

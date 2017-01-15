function x = apply_Tfull(obj,x0,k)
% Do x = (A - k^2*B)*x0 and enforce Dirichlet BCs at outer edge.
% Inputs:
%  - obj: A PML instance
%  - x0 : vector we're applying obj.T(k) to quickly
%  - k  : energy parameter, i.e. T(k) = A - k^2*B

x = obj.apply_op_no_bcs(x0,k^2,obj.VvaluesVec);

% Dirichlet boundary condition row action
x(obj.BCrows) = x0(obj.BCrows);

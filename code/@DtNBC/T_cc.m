function Tk = T_cc(obj,k)
% computes matrix-valued function T whose eigs are the resonances
    Tk = obj.A_cc - k^2*obj.B_cc + obj.C_cc(k);
end
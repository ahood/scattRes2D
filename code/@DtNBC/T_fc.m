function Tk = T_fc(obj,k)
% computes matrix-valued function T whose eigs are the resonances
    Tk = obj.A_fc - k^2*obj.B_fc + obj.C_fc(k);
end
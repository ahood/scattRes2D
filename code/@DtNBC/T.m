function Tk = T(obj,k)
% computes matrix-valued function T whose eigs are the resonances
    Tk = obj.A - k^2*obj.B + obj.C(k);
end
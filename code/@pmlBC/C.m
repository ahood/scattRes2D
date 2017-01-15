function Ck = C(obj,k)
    Ck = obj.pieces.update(obj.A - k^2*obj.B);
end
function Tk = Tfull(obj,k)

    % the rational approximation is expecting an energy, not k
    z = k^2;
    Tk = obj.A - z*obj.B;
end
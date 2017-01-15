function dfoverf = dLogTz(obj,z,br)
% implements df/f where f = det T, and
% T(z) = obj.T(sqrt(z)), with the br-th branch
% of sqrt
    if br == 1, kz = sqrt(z); else kz = -sqrt(z); end
    dkdz = kz^(-1)/2;
    Tz = obj.T(kz);
    dTdk = obj.dT(kz);
    dTdz = dTdk*dkdz;
    dfoverf = trace( Tz\dTdz );
end
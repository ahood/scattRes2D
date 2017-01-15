function d_dfoverf = d2LogTz(obj,z,br)
% derivative of dLogTz
    if br == 1, kz = sqrt(z); else, kz = -sqrt(z); end
    dkdz = kz^(-1)/2;
    d2kdz = -kz^(-2)*dkdz/2;
    Tz = obj.T(kz);
    dTdk = obj.dT(kz);
    dTdz = dTdk*dkdz;
    d2Tdk = obj.d2T(kz);
    d2Tdz = d2Tdk*dkdz^2 + dTdk*d2kdz;
    d_dfoverf = trace( -(Tz\dTdz)^2 + Tz\d2Tdz );
end
function ratApproxPoles = get_poles_in_range(n,region,R,br)
    addpath(fileparts(pwd));

    n = abs(n);
    filename = sprintf('../hankel_zeros/n_%03d.dat',n);
    tmp = csvread(filename);
    ws = tmp(:,1) + 1i*tmp(:,2); % poles of hankel functions
    % pick out the ones in range of given branch of sqrt function
    if br == 1 % principal branch, range RHP
        ws = ws( real(ws) > 0 | (real(ws) == 0 & imag(ws) > 0) );
    else
        ws = ws( real(ws) < 0 | (real(ws) == 0 & imag(ws) < 0) );
    end
    ks = ws/R;
    zs = ks.^2; % poles of DtN map in z-plane
    zs = zs( region.contains(zs) ); % keep the ones in given region
    rounded = 1e-10*round(1e10*zs);
    ratApproxPoles = unique(rounded);
end
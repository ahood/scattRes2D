function save_hankel_zeros

% close all

nmax = 500;
smax = 500; % number zeros to compute for each n

% figure, hold on
for n = 0:nmax
    zns         = hankel_zeros_large(n,smax,[],[]);
    [azns,bzns] = hankel_zeros_small(n,smax,[],[]);
    filename = sprintf('n_%03d.dat', n);
    f = fopen(filename,'w');
    for z = [zns, azns, bzns] % didn't bother to sort
        fprintf(f, '%+14.10e, %+14.10e\n', real(z), imag(z) );
    end
    fclose(f);
    
%     plot(real( zns),imag( zns),'*k');
%     plot(real(azns),imag(azns),'*b');
%     plot(real(bzns),imag(bzns),'*r');
end

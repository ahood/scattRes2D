function show_truescattn_sizes(maxn,k,Vs,s)

truescattFourier = zeros(2*maxn+1);
fourier_ns = -maxn:maxn;
for n = fourier_ns
    truescattFourier(n+maxn+1) = norm(compute_truescattn(n,k,Vs,s));
end
figure, plot(fourier_ns,log10(truescattFourier));
title('log10 norm of fourier coefficients');

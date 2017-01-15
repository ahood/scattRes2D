function test_scattResComp2d_plotting(verbose)

close all

if nargin == 0
    verbose = 0;
end

% set up a generic scattering problem 
Nt = 30; Nr = 50; R = 2;
s = scattResComp2d_sparse(Nt,Nr,R);

% pick a function and express in rectangular, polar, and complex
% coordinates
fRect = @(x,y) (1-y).*exp(1i*5*x);
% fRect = @(x,y) exp(5i*x);
fPolar = @(r,t) fRect(r.*cos(t),r.*sin(t));
fComplex = @(z) fRect(real(z),imag(z));

% compare (numerically) the functions evaluated on the mesh
part = @real;
uRectR    = s.plotFun(fRect   ,'real(uRect)'   ,part,'rect'   ); close
uPolarR   = s.plotFun(fPolar  ,'real(uPolar)'  ,part,'polar'  ); close 
uComplexR = s.plotFun(fComplex,'real(uComplex)',part,'complex'); close
abserr = norm( uRectR - uPolarR ); relerr = abserr/norm( uRectR );
if verbose | abserr > 1e-16
    fprintf('Error between real parts of uRect and uPolar: %4.2e, %4.2e\n', abserr,relerr);
end
abserr = norm( uRectR - uComplexR ); relerr = abserr/norm( uRectR );
if verbose | abserr > 1e-16
    fprintf('Error between real parts of uRect and uComplex: %4.2e, %4.2e\n', abserr,relerr);
end

part = @imag;
uRectI    = s.plotFun(fRect   ,'imag(uRect)'   ,part,'rect'   ); close
uPolarI   = s.plotFun(fPolar  ,'imag(uPolar)'  ,part,'polar'  ); close
uComplexI = s.plotFun(fComplex,'imag(uComplex)',part,'complex'); close
abserr = norm( uRectI - uPolarI ); relerr = abserr/norm( uRectI );
if verbose | abserr > 1e-16
    fprintf('Error between imag parts of uRect and uPolar: %4.2e, %4.2e\n', abserr,relerr);
end
abserr = norm( uRectI - uComplexI ); relerr = abserr/norm( uRectI );
if verbose | abserr > 1e-16
    fprintf('Error between imag parts of uRect and uComplex: %4.2e, %4.2e\n', abserr,relerr);
end

% Make sure valuesVec same for each kind of coordinates
valuesVecRect    = s.valuesVecFromFun(fRect   ,'rect'   );
valuesVecPolar   = s.valuesVecFromFun(fPolar  ,'polar'  );
valuesVecComplex = s.valuesVecFromFun(fComplex,'complex');
abserr = norm( valuesVecRect - valuesVecPolar ); relerr = abserr/norm( valuesVecRect );
if verbose | abserr > 1e-16
    fprintf('Error between valuesVecRect and valuesVecPolar: %4.2e, %4.2e\n', abserr,relerr);
end
abserr = norm( valuesVecRect - valuesVecComplex ); relerr = abserr/norm( valuesVecRect );
if verbose | abserr > 1e-16
    fprintf('Error between valuesVecRect and valuesVecComplex: %4.2e, %4.2e\n', abserr,relerr);
end

% Check that valuesVec and fourierVec have same plots
fourierVecRect = s.fourierVecFromFun(fRect,'rect');

part = @real;
uFourierVecRectR = s.plotFourierVec(fourierVecRect,'',part); close
uValuesVecRectR  = s.plotValuesVec(valuesVecRect,'',part);   close
abserr = norm( uFourierVecRectR - uValuesVecRectR ); relerr = abserr/norm( uValuesVecRectR );
if verbose | abserr > 1e-12
    fprintf('Error in real parts of plotted valuesVec and fourierVec: %4.2e, %4.2e\n', abserr,relerr);
end

part = @imag;
uFourierVecRectI = s.plotFourierVec(fourierVecRect,'',part); close
uValuesVecRectI  = s.plotValuesVec(valuesVecRect,'',part);   close
abserr = norm( uFourierVecRectI - uValuesVecRectI ); relerr = abserr/norm( uValuesVecRectI );
if verbose | abserr > 1e-12
    fprintf('Error in imag parts of plotted valuesVec and fourierVec: %4.2e, %4.2e\n', abserr,relerr);
end

% check that plot from function and its valuesVec are same
abserr = norm( uValuesVecRectR - uRectR ); relerr = abserr/norm( uRectR );
if verbose | abserr > 1e-16
    fprintf('Error in real parts of plotted function and valuesVec: %4.2e, %4.2e\n', abserr,relerr);
end
abserr = norm( uValuesVecRectI - uRectI ); relerr = abserr/norm( uRectI );
if verbose | abserr > 1e-16
    fprintf('Error in imag parts of plotted function and valuesVec: %4.2e, %4.2e\n', abserr,relerr);
end

% Check that RHS same for each kind of coordinates
VRect   = @(x,y) x.^2.*(y - x.*y.^2);
dir     = dirBC(s.Nt,s.Nrs,{VRect},'rect',s.Rs); % cleaned up RHSfromFun needs BC instance
RHSRect = dir.RHSfromFun(fRect);

VPolar   = @(r,t) VRect(r.*cos(t),r.*sin(t));
dir      = dirBC(s.Nt,s.Nrs,{VPolar},'polar',s.Rs);
RHSPolar = dir.RHSfromFun(fPolar);

VComplex   = @(z) VRect(real(z),imag(z));
dir        = dirBC(s.Nt,s.Nrs,{VComplex},'complex',s.Rs);
RHSComplex = dir.RHSfromFun(fComplex);

abserr = norm( RHSRect - RHSPolar ); relerr = abserr/norm( RHSRect );
if verbose | abserr > 1e-16
    fprintf('Error between RHSRect and RHSPolar: %4.2e, %4.2e\n', abserr,relerr);
end
abserr = norm( RHSRect - RHSComplex ); relerr = abserr/norm( RHSRect );
if verbose | abserr > 1e-16
    fprintf('Error between RHSRect and RHSComplex: %4.2e, %4.2e\n', abserr,relerr);
end

% Check that fourierVec is what it should be if function has finite Fourier
% expansion
Nt = 4; Nr = 50; R = 2;
s = scattResComp2d_sparse(Nt,Nr,R);
fm1 = @(r) r.^2;
f0  = @(r) cos(r);
fp1 = @(r) 1 + r;
f = @(r,t) fm1(r).*exp(1i*(-1)*t) + f0(r) + fp1(r).*exp(1i*1*t);
fhat = @(r) [fm1(r), f0(r), fp1(r), 0*r];

fourierVec = s.fourierVecFromFun(f,'polar');
exactFourierVec = [fm1(s.r); f0(s.r); fp1(s.r); 0*s.r];

abserr = norm( fourierVec - exactFourierVec ); relerr = abserr/norm( fourierVec );
if verbose | abserr > 1e-14
    fprintf('Error between fourierVec and what it should be (finite exp): %4.2e, %4.2e\n', abserr,relerr);
end

% Check that things working properly if f = exp(1i*k*x), which has
% Jacobi-Anger expansion
Nt = 60; Nr = 50; R = 2;
genericScatt = scattResComp2d(Nt,Nr,R);
k = 5;
f = @(r,t) exp(1i*k*r.*cos(t));
fc = @(n,z) 1i.^n.*besselj(n,z); % n-th fourier coeff of f
t = genericScatt.theta(2); % try the first ray off the real axis
fapprox = 0*genericScatt.r;
for n = genericScatt.Ns
    fapprox = fapprox + fc(n,k*genericScatt.r)*exp(1i*n*t);
end
abserr = norm( fapprox - f(genericScatt.r,t) ); relerr = abserr/norm( fapprox );
if verbose | abserr > 1e-11
    fprintf('Error between true e^{ikx} and approx reconstructed from fourier exp: %4.2e, %4.2e\n',...
    abserr,relerr); 
end

longr = repmat(genericScatt.r,1,genericScatt.Nt); longr = longr(:);
% longn = repmat(-Ns:Ns,genericScatt.Nr,1);         longn = longn(:);
longn = repmat(genericScatt.Ns,genericScatt.Nr,1); longn = longn(:);
exactFourierVec = fc(longn,k*longr);
fourierVec = genericScatt.fourierVecFromFun(f,'polar');
abserr = norm( exactFourierVec - fourierVec ); relerr = abserr/norm( exactFourierVec );
if verbose | abserr > 1e-11
    fprintf('Error between automatically created fourierVec and exact one: %4.2e, %4.2e\n', ...
    abserr,relerr);
end

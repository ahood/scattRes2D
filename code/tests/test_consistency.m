function test_consistency(verbose)

addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end

%% Default value of incfun
% The user can pass incfun or a wave number k, and in the latter case
% incfun is defined to be exp(1i*k*x). Check each gives the same result.
Nt = 10;
Nrs = 30;
Rs = 3;
Vs = {@(r,t) cos(r).*sin(t)};
coords = 'polar';
k1 = rand;
incfun = @(r,t) exp(1i*k1*r.*cos(t));

% check dirBC
dir = dirBC(Nt,Nrs,Vs,coords,Rs);
RHS_dir = dir.RHSfromFun(incfun);
[abserr,relerr] = get_err(RHS_dir, dir.RHSfromFun(k1));
if verbose || abserr > 1e-10
    fprintf('Error in dirBC RHS created with default incfun: %4.2e, %4.2e\n', abserr, relerr);
end

% check pmlBC
Nr_out = 10;
Vs_pml = [Vs, {@(r,t) 0*r}];
Rout = Rs(end) + 1;
a = 0.4;
l = @(t) a*(t-Rout).^3;
dl = @(t) a*3*(t-Rout).^2;
d2l = @(t) a*3*2*(t-Rout);
pml = pmlBC(Nt,Nrs,Nr_out,Vs_pml,coords,Rs,Rout,l,dl,d2l);
RHS_pml = pml.RHSfromFun(incfun);
[abserr,relerr] = get_err(RHS_pml, pml.RHSfromFun(k1));
if verbose || abserr > 1e-10
    fprintf('Error in pmlBC RHS created with default incfun: %4.2e, %4.2e\n', abserr, relerr);
end

% check pmlBC full version (finds scatt wave on PML region as well)
RHS_pml_full = pml.RHSfromFun_full(incfun);
[abserr,relerr] = get_err(RHS_pml_full, pml.RHSfromFun_full(k1));
if verbose || abserr > 1e-10
    fprintf('Error in pmlBC full RHS created with default incfun: %4.2e, %4.2e\n', abserr, relerr);
end

% check DtNBC
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
RHS_dtn = dtn.RHSfromFun(incfun);
[abserr,relerr] = get_err(RHS_dtn, dtn.RHSfromFun(k1));
if verbose || abserr > 1e-10
    fprintf('Error in DtNBC RHS created with default incfun: %4.2e, %4.2e\n', abserr, relerr);
end

% make sure RHS same in dirBC, pmlBC, and DtNBC cases
[abserr,relerr] = get_err(RHS_dtn,RHS_dir);
if verbose || abserr > 1e-10
    fprintf('Error between DtNBC and dirBC RHS: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = get_err(RHS_dtn,RHS_pml);
if verbose || abserr > 1e-10
    fprintf('Error between DtNBC and pmlBC RHS: %4.2e, %4.2e\n', abserr, relerr);
end

% make sure the part of RHS not in PML region matches with other RHS's
RHS_pml_full2 = reshape(RHS_pml_full,pml.Nr,pml.Nt);
RHS_pml_full2(end-pml.Nrs(end)+1:end,:) = [];
RHS_pml_full2 = RHS_pml_full2(:);
[abserr,relerr] = get_err(RHS_pml,RHS_pml_full2);
if verbose || abserr > 1e-10
    fprintf('Error between pmlBC RHS and full RHS on non-PML region: %4.2e, %4.2e\n', abserr, relerr);
end

%% Converting between representations
% Testing
% 1. valuesVecFromFourierVec and fourierVecFromValuesVec
% 2. valuesVecFromFun and fourierVecFromFun
% 3. valuesVecFromFunCellArray and fourierVecFromFunCellArray
% 4. chebVecFromFourierVec and fourierVecFromChebVec

% set up instance
Nt = 4; Nrs = [20,25]; Rs = [1.1,2];
s = scattResComp2d(Nt,Nrs,Rs);

% pick a function for evaluation, defined first in term of its fourier
% coefficient functions
fm1 = @(r) r.^2 + 5;
f0  = @(r) r;
fp1 = @(r) 1./r - 1;
fp2 = @(r) cos(r);

% for special case where more fourier coeffs than points in theta mesh
fm2 = @(r) 0*r;
fp3 = @(r) 0*r; 

% function
f = @(r,t) fm1(r).*exp(-1i*t) + f0(r) + fp1(r).*exp(1i*t) + fp2(r).*exp(2i*t);

% unnecessarily piecewise definition
fs = {f,f};

% exact valuesVec and fourierVec
exactValuesVec = f(s.rs,s.ts);
exactFourierVec = [fm1(s.r); f0(s.r); fp1(s.r); fp2(s.r)];
longExactFourierVec = [fm2(s.r); exactFourierVec; fp3(s.r)];

% different coordinates
frect = @(x,y) f(sqrt(x.^2 + y.^2), atan2(y,x));
fcomp = @(z)   f(abs(z), angle(z));

%%
% 1. valuesVecFromFourierVec and fourierVecFromValuesVec
[abserr,relerr] = get_err(exactValuesVec, s.valuesVecFromFourierVec(exactFourierVec));
if verbose || abserr > 1e-13
    fprintf('Error for valuesVecFromFourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(exactFourierVec, s.fourierVecFromValuesVec(exactValuesVec));
if verbose || abserr > 1e-14
    fprintf('Error for fourierVecFromValuesVec: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(exactValuesVec, s.valuesVecFromFourierVec(longExactFourierVec));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFourierVec (more coeffs than thetas): %4.2e, %4.2e\n', abserr, relerr);
end

%%
% 2. valuesVecFromFun and fourierVecFromFun
[abserr,relerr] = get_err(exactValuesVec, s.valuesVecFromFun(f,'polar'));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFun (polar):   %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(exactValuesVec, s.valuesVecFromFun(frect,'rect'));
if verbose || abserr > 1e-13
    fprintf('Error for valuesVecFromFun (rect):    %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(exactValuesVec, s.valuesVecFromFun(fcomp,'complex'));
if verbose || abserr > 1e-13
    fprintf('Error for valuesVecFromFun (complex): %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(exactFourierVec, s.fourierVecFromFun(f,'polar'));
if verbose || abserr > 1e-14
    fprintf('Error for fourierVecFromFun (polar):  %4.2e, %4.2e\n', abserr, relerr);
end

%%
% 3. valuesVecFromFunCellArray and fourierVecFromFunCellArray
[abserr,relerr] = get_err(exactValuesVec, s.valuesVecFromFunCellArray(fs,'polar'));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFunCellArray (polar): %4.2e, %4.2e\n', abserr, relerr); 
end

[abserr,relerr] = get_err(exactFourierVec, s.fourierVecFromFunCellArray(fs,'polar'));
if verbose || abserr > 1e-14
    fprintf('Error for fourierVecFromFunCellArray (polar): %4.2e, %4.2e\n', abserr, relerr);
end

%%
% 4. chebVecFromFourierVec and fourierVecFromChebVec
% Pick Cheb coeffs and fourierVec, then find chebVec

% set up instance
Nt = 20; Nrs = [20,25]; Rs = [1.1,2];
s = scattResComp2d(Nt,Nrs,Rs);

% make up some chebyshev coefficients for piecewise-defined fourier coeffs
chebRect1 = zeros(Nrs(1),Nt);
chebRect2 = zeros(Nrs(2),Nt);
% even columns for even functions since Nt/2 even
% use f(x) = x^2 as underlying function, expansion 0.5*T0 + 0.5*T2
chebRect1([1,2],2:2:end) = 1/2; % only take even poly coeffs
chebRect2([1,3],2:2:end) = 1/2;
% odd columns for odd functions
% use f(x) = x as underlying function, expansion T1
chebRect1(1,1:2:end) = 1; % only take odd poly coeffs
chebRect2(2,1:2:end) = 1;
chebRect = [chebRect1; chebRect2];

% manually create the values of fourier coeffs
% will be expansions of above functions on meshes of [-1,1] ([0,1] for first)
[~,x1] = cheb(2*Nrs(1)-1); x1 = x1(Nrs(1)+1:end);
[~,x2] = cheb(  Nrs(2)-1);
x = [x1; x2]; % evaluate on this, NOT on s.r
v_even = x.^2;
v_odd  = x;
fourierRect = zeros(s.Nr,s.Nt);
fourierRect(:,2:2:end) = repmat(v_even, 1, Nt/2);
fourierRect(:,1:2:end) = repmat(v_odd,  1, Nt/2);
fourierVec = fourierRect(:);

% make chebVec from fourierVec and compare
chebVec = s.chebVecFromFourierVec(fourierVec);

[abserr,relerr] = get_err(chebRect(:), chebVec);
if verbose || abserr > 1e-13
    fprintf('Error for chebVecFromFourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

% check that I kron Winv maps cheb coeffs to values
chebRect = reshape(chebVec, s.Nr, s.Nt);
fourierRect = 0*chebRect;
if mod(s.Nt/2,2) == 0 % last column even
    cols_even = 2:2:s.Nt;
    cols_odd  = 1:2:s.Nt;
else % last column odd
    cols_even = 1:2:s.Nt;
    cols_odd  = 2:2:s.Nt;
end
fourierRect(:,cols_even) = s.Winv_even*chebRect(:,cols_even);
fourierRect(:,cols_odd ) = s.Winv_odd *chebRect(:,cols_odd );

[abserr,relerr] = get_err(fourierVec, fourierRect(:));
if verbose || abserr > 1e-12
    fprintf('Error in Winv action: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% check that Winv_even/odd is the inverse of radialVals2chebCoeffs

% radial vals of even fcn
f = @(x) x.^2;
v = f(s.r);

[abserr,relerr] = get_err(v,s.Winv_even*s.radialVals2chebCoeffs(v,'even'));
if verbose || abserr > 1e-12
    fprintf('Error in Winv_even*radialVals2chebCoeffs = I for f(x) = x^2: %4.2e, %4.2e\n', abserr, relerr);
end

% radial vals of even fcn of higher degree
f = @(x) cos(x);
v = f(s.r);

[abserr,relerr] = get_err(v,s.Winv_even*s.radialVals2chebCoeffs(v,'even'));
if verbose || abserr > 1e-13
    fprintf('Error in Winv_even*radialVals2chebCoeffs = I for f(x) = cos(x): %4.2e, %4.2e\n', abserr, relerr);
end

% radial vals of odd fcn
f = @(x) x;
v = f(s.r);

[abserr,relerr] = get_err(v,s.Winv_odd*s.radialVals2chebCoeffs(v,'odd'));
if verbose || abserr > 1e-13
    fprintf('Error in Winv_odd*radialVals2chebCoeffs = I for f(x) = x: %4.2e, %4.2e\n', abserr, relerr);
end

% radial vals of odd fcn of higher degree
f = @(x) sin(x);
v = f(s.r);

[abserr,relerr] = get_err(v,s.Winv_odd*s.radialVals2chebCoeffs(v,'odd'));
if verbose || abserr > 1e-13
    fprintf('Error in Winv_odd*radialVals2chebCoeffs = I for f(x) = sin(x): %4.2e, %4.2e\n', abserr, relerr);
end

%% T(k) in the Chebyshev basis

Nt = 20; Nrs = [20,30];

%%
% with Dirichlet BCs
dir = dirBC(Nt,Nrs,Vs,coords,Rs);

scattValuesVec_dir  = dir.solve(k1);
scattFourierVec_dir = dir.solve_fc(k1);
scattChebVec_dir    = dir.solve_cc(k1);

[abserr,relerr] = get_err(scattFourierVec_dir, dir.fourierVecFromValuesVec(scattValuesVec_dir));
if verbose || abserr > 1e-12
    fprintf('Error in scattValuesVec vs scattFourierVec for dir BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(scattChebVec_dir,    dir.chebVecFromValuesVec(scattValuesVec_dir));
if verbose || abserr > 1e-14
    fprintf('Error in scattValuesVec vs scattChebVec for dir BCs: %4.2e, %4.2e\n', abserr, relerr);
    fprintf('  Can''t seem to get better than relerr 1e-4...\n');
end

%%
% with PML BCs
Nr_out = 10;
Vs_pml = [Vs, {@(r,t) 0*r}];
Rout = Rs(end) + 1;
a = 0.4;
l = @(t) a*(t-Rout).^3;
dl = @(t) a*3*(t-Rout).^2;
d2l = @(t) a*3*2*(t-Rout);
pml = pmlBC(Nt,Nrs,Nr_out,Vs_pml,coords,Rs,Rout,l,dl,d2l);

scattValuesVec_pml  = pml.solve(k1);
scattFourierVec_pml = pml.solve_fc(k1);
scattChebVec_pml    = pml.solve_cc(k1);

[abserr,relerr] = get_err(scattFourierVec_pml, dir.fourierVecFromValuesVec(scattValuesVec_pml));
if verbose || abserr > 1e-10
    fprintf('Error in scattValuesVec vs scattFourierVec for pml BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(scattChebVec_pml,    dir.chebVecFromValuesVec(scattValuesVec_pml));
if verbose || abserr > 1e-14
    fprintf('Error in scattValuesVec vs scattChebVec for pml BCs: %4.2e, %4.2e\n', abserr, relerr);
    fprintf('  Can''t seem to get better than relerr 1e-4...\n');
end

scattValuesVec_pml_full  = pml.solve_full(k1);
scattFourierVec_pml_full = pml.solve_full_fc(k1);
scattChebVec_pml_full    = pml.solve_full_cc(k1);

[abserr,relerr] = get_err(scattFourierVec_pml_full, pml.fourierVecFromValuesVec(scattValuesVec_pml_full));
if verbose || abserr > 1e-10
    fprintf('Error in scattValuesVec vs scattFourierVec for pml with dir BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(scattChebVec_pml_full,    pml.chebVecFromValuesVec(scattValuesVec_pml_full));
if verbose || abserr > 1e-14
    fprintf('Error in scattValuesVec vs scattChebVec for pml with dir BCs: %4.2e, %4.2e\n', abserr, relerr);
    fprintf('  Can''t seem to get better than relerr 1e-4...\n');
end

[abserr,relerr] = get_err(scattValuesVec_pml, scattValuesVec_pml_full(pml.pieces.I1));
if verbose || abserr > 1e-11
    fprintf('Error between pml solutions using extended and usual mesh: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% with DtN BCs
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);

scattValuesVec_dtn  = dtn.solve(k1);
scattFourierVec_dtn = dtn.solve_fc(k1);
scattChebVec_dtn    = dtn.solve_cc(k1);

[abserr,relerr] = get_err(scattFourierVec_dtn, dtn.fourierVecFromValuesVec(scattValuesVec_dtn));
if verbose || abserr > 1e-12
    fprintf('Error in scattValuesVec vs scattFourierVec for DtN BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(scattChebVec_dtn,    dtn.chebVecFromValuesVec(scattValuesVec_dtn));
if verbose || abserr > 1e-14
    fprintf('Error in scattValuesVec vs scattChebVec for DtN BCs: %4.2e, %4.2e\n', abserr, relerr);
    fprintf('  Can''t seem to get better than relerr 1e-4...\n');
end

%% Internal consistency of dtnBC

% choose piecewise constant, axisymmetric 2D potential
Rin = 1; Rout = 2; V0_in = 0; V0_out = 50;
coords = 'polar';
Vs = {@(r,t) 0*r + V0_in; @(r,t) 0*r + V0_out};

% choose discretization parameters and PML parameters
Nt = 16; Nrs = [20,20]; Rs = [Rin,Rout];
RPML = 6; Nr_out = 40;
a = 0.4;
l = @(t) a*(t-Rout).^3;
dl = @(t) a*3*(t-Rout).^2;
d2l = @(t) a*3*2*(t-Rout);
pmlVs = Vs; pmlVs{end+1} = @(r,t) 0*r;

% make dtn problem
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);

% test derivatives of DtN map
k0 = 5*rand + 5i*rand;
dk0 = 1e-6*rand + 1i*1e-6*rand;
n = randi([-10,10]);
Rtest = 10*rand + 0.5;

f   = @(k) DtNBC.DtNcoeffs(n,k,Rtest);
df  = @(k) DtNBC.dDtNcoeffs(n,k,Rtest);
d2f = @(k) DtNBC.d2DtNcoeffs(n,k,Rtest);

[abserr,relerr] = get_deriv_err(f,df,k0,dk0);
if verbose || abserr > 1e-8
    fprintf('Consistency error between f and df: %4.2e, %4.2e\n', ...
        abserr,relerr);
end

[abserr,relerr] = get_deriv_err(df,d2f,k0,dk0);
if verbose || abserr > 1e-7
    fprintf('Consistency error between df and d2f: %4.2e, %4.2e\n', ...
        abserr,relerr);
end

% test derivatives of T
[abserr,relerr] = get_deriv_err(@dtn.T,@dtn.dT,k0,dk0);
if verbose || abserr > 1e-5
    fprintf('Consistency error between T(k) and dT(k): %4.2e, %4.2e\n', ...
        abserr,relerr);
end

[abserr,relerr] = get_deriv_err(@dtn.dT,@dtn.d2T,k0,dk0);
if verbose || abserr > 1e-7
    fprintf('Consistency error between dT(k) and d2T(k): %4.2e, %4.2e\n', ...
        abserr,relerr);
end

% check df/f where f = det T
    function kz = k(z)
        kz = sqrt(z);
    end
    function dkdz = dk(z)
        dkdz = k(z)^(-1)/2;
    end
    function d2kdz = d2k(z)
        d2kdz = -k(z).^(-2).*dk(z)/2;
    end  
    function Tz = T(z)
        kz = k(z);
        Tz = dtn.T(kz);
    end
    function dTdz = dT(z)
        kz = k(z);
        dTdk = dtn.dT(kz);
        dTdz = dTdk*dk(z);
    end
    function d2Tdz = d2T(z)
        kz = k(z);
        dTdk = dtn.dT(kz); d2Tdk = dtn.d2T(kz);
        d2Tdz = d2Tdk*dk(z)^2 + dTdk*d2k(z);
    end

k0 = 10*rand + 10i*rand;
dk0 = 1e-6*rand + 1i*1e-6*rand;
z0 = 10*rand + 10i*rand;
dz = 1e-6*rand + 1i*1e-6*rand;

[abserr,relerr] = get_deriv_err(@k,@dk,z0,dz);
if verbose || abserr > 1e-9
    fprintf('Consistency between k(z) and dk(z): %4.2e, %4.2e\n', ...
        abserr,relerr);
end
[abserr,relerr] = get_deriv_err(@dk,@d2k,z0,dz);
if verbose || abserr > 1e-9
    fprintf('Consistency between dk(z) and d2k(z): %4.2e, %4.2e\n', ...
        abserr,relerr);
end

[abserr,relerr] = get_deriv_err(@T,@dT,z0,dz);
if verbose || abserr > 1e-5
    fprintf('Consistency between T(z) and dT(z): %4.2e, %4.2e\n', ...
        abserr,relerr);
end
[abserr,relerr] = get_deriv_err(@dT,@d2T,z0,dz);
if verbose || abserr > 1e-8
    fprintf('Consistency between dT(z) and d2T(z): %4.2e, %4.2e\n', ...
        abserr,relerr);
end

dfoverf = @(z) trace( T(z)\dT(z) );
d_dfoverf = @(z) trace( -(T(z)\dT(z))^2 + T(z)\d2T(z) );

[abserr,relerr] = get_deriv_err(dfoverf,d_dfoverf,z0,dz);
if verbose || abserr > 1e-8
    fprintf('Consistency between dfoverf and deriv: %4.2e, %4.2e\n', ...
        abserr,relerr);
end

[abserr,relerr] = get_err(dfoverf(z0),dtn.dLogTz(z0,1));
if verbose || abserr > 1e-16
    fprintf('Error in DtNBC.dLogTz: %4.2e, %4.2e\n', abserr,relerr);
end
[abserr,relerr] = get_err(d_dfoverf(z0),dtn.d2LogTz(z0,1));
if verbose || abserr > 1e-16
    fprintf('Error in DtNBC.d2LogTz: %4.2e, %4.2e\n', abserr,relerr);
end

%% pml.T and dtn.T

pml = pmlBC(Nt,Nrs,Nrs(end),pmlVs,coords,Rs,RPML,l,dl,d2l);
k0 = 10*rand;
[abserr,relerr] = get_err(dtn.T(k0),pml.T(k0));
if verbose || relerr > 1e-2
    fprintf('Consistency error between dtn.T and pml.T: %4.2e, %4.2e\n', ...
        abserr,relerr);
end

[abserr,relerr] = get_deriv_err(@pml.T,@pml.dT,k0,dk0);
if verbose || relerr > 1e-7
    fprintf('Consistency error between pml.T and pml.dT: %4.2e, %4.2e\n', ...
        abserr,relerr);
end

end

function [abserr,relerr] = get_deriv_err(f,df,x,dx)
    dfapprox = (f(x+dx)-f(x-dx))/2/dx;
    [abserr,relerr] = get_err(df(x),dfapprox);
end

function [abserr,relerr] = get_err(x_true,x_approx)
    abserr = norm(full(x_true-x_approx));
    relerr = abserr/norm(full(x_true));
end
function test_compute_truescatt(verbose)
% For several piecewise constant potentials, the true solution to the 
% scattering problem is computed and compared to the solution computed by
% using DtN boundary conditions.

close all
addpath(fileparts(pwd));
disp('Temporarily turning warnings off');
warning off % inverting nearly singular matrices to solve

if nargin == 0
    verbose = 0;
end

% set down parameters common to all tests
k = 2; incfun = @(r,t) exp(1i*k*r.*cos(t)); % incident wave
abstol = 1e-8; 
reltol = 1e-8;

%% Test 1: two subregions, Rs = [1,2], V0s = [0,50]

% pick axisymmetric potential, piecewise constant on concentric annuli
Rs = [1,2]; V0s = [0,50]; 
N = length(V0s); coords = 'polar';
Vs = cell(1,N); for j = 1:N, Vs{j} = @(r,t) 0*r + V0s(j); end

% choose discretization parameters, make DtNBC object, solve
Nt = 30; Nrs = 0*V0s + 20;
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtnscatt = dtn.solve(k);

% figure out how many fourier coeffs we need
maxn = 24;
truescatt = compute_truescatt(k,V0s,dtn,maxn);
truescatt2 = compute_truescatt(k,V0s,dtn,2*maxn);
[abserr,relerr] = L2err(dtn,truescatt,dtn,truescatt2,[dtn.theta, 2*pi]);

% collect maximum absolute and relative errors
[abserr,relerr] = L2err(dtn,truescatt,dtn,dtnscatt,[dtn.theta, 2*pi]);
if verbose | abserr > 1e-10
    fprintf('L2 error for Rs = [1,2], V0s = [0,50]: %4.2e, %4.2e\n', abserr, relerr);
end

%% Test 2: two subregions, Rs = [1,2], V0s = [20,50]

% pick axisymmetric potential, piecewise constant on concentric annuli
Rs = [1,2]; V0s = [20,50]; 
N = length(V0s); coords = 'polar';
Vs = cell(1,N); for j = 1:N, Vs{j} = @(r,t) 0*r + V0s(j); end

% choose discretization parameters, make DtNBC object, solve
Nt = 30; Nrs = 0*V0s + 20;
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtnscatt = dtn.solve(k);

% figure out how many fourier coeffs we need
maxn = 24;
truescatt = compute_truescatt(k,V0s,dtn,maxn);
truescatt2 = compute_truescatt(k,V0s,dtn,2*maxn);
[abserr,relerr] = L2err(dtn,truescatt,dtn,truescatt2,[dtn.theta, 2*pi]);

% collect maximum absolute and relative errors
[abserr,relerr] = L2err(dtn,truescatt,dtn,dtnscatt,[dtn.theta, 2*pi]);
if verbose | abserr > 1e-10
    fprintf('L2 error for Rs = [1,2], V0s = [0,50]: %4.2e, %4.2e\n', abserr, relerr);
end

%% Test 3: two subregions, Rs = [1.5,2.75], V0s = [20,50]

% pick axisymmetric potential, piecewise constant on concentric annuli
Rs = [1,2]; V0s = [20,50]; 
N = length(V0s); coords = 'polar';
Vs = cell(1,N); for j = 1:N, Vs{j} = @(r,t) 0*r + V0s(j); end

% choose discretization parameters, make DtNBC object, solve
Nt = 30; Nrs = 0*V0s + 20;
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtnscatt = dtn.solve(k);

% figure out how many fourier coeffs we need
maxn = 24;
truescatt = compute_truescatt(k,V0s,dtn,maxn);
truescatt2 = compute_truescatt(k,V0s,dtn,2*maxn);
[abserr,relerr] = L2err(dtn,truescatt,dtn,truescatt2,[dtn.theta, 2*pi]);

% collect maximum absolute and relative errors
[abserr,relerr] = L2err(dtn,truescatt,dtn,dtnscatt,[dtn.theta, 2*pi]);
if verbose | abserr > 1e-10
    fprintf('L2 error for Rs = [1,2], V0s = [0,50]: %4.2e, %4.2e\n', abserr, relerr);
end

%% Test 4: three subregions, Rs = [1,2,3], V0s = [20,50,30]

% pick axisymmetric potential, piecewise constant on concentric annuli
Rs = [1,2,3]; V0s = [20,50,30]; 
N = length(V0s); coords = 'polar';
Vs = cell(1,N); for j = 1:N, Vs{j} = @(r,t) 0*r + V0s(j); end

% choose discretization parameters, make DtNBC object, solve
Nt = 40; Nrs = 0*V0s + 20;
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtnscatt = dtn.solve(k);

% figure out how many fourier coeffs we need
maxn = 25;
truescatt = compute_truescatt(k,V0s,dtn,maxn);
truescatt2 = compute_truescatt(k,V0s,dtn,2*maxn);
[abserr,relerr] = L2err(dtn,truescatt,dtn,truescatt2,[dtn.theta, 2*pi]);

% collect maximum absolute and relative errors
[abserr,relerr] = L2err(dtn,truescatt,dtn,dtnscatt,[dtn.theta, 2*pi]);
if verbose | abserr > 1e-8
    fprintf('L2 error for Rs = [1,2], V0s = [0,50]: %4.2e, %4.2e\n', abserr, relerr);
end

%% Test 5: just one region, Rs = 1.5, V0s = 10

% pick axisymmetric potential, piecewise constant on concentric annuli
Rs = 1.5; V0s = 10; 
N = length(V0s); coords = 'polar';
Vs = cell(1,N); for j = 1:N, Vs{j} = @(r,t) 0*r + V0s(j); end

% choose discretization parameters, make DtNBC object, solve
Nt = 30; Nrs = 0*V0s + 20;
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtnscatt = dtn.solve(k);

% figure out how many fourier coeffs we need
maxn = 25;
truescatt = compute_truescatt(k,V0s,dtn,maxn);
truescatt2 = compute_truescatt(k,V0s,dtn,2*maxn);
[abserr,relerr] = L2err(dtn,truescatt,dtn,truescatt2,[dtn.theta, 2*pi]);

% collect maximum absolute and relative errors
[abserr,relerr] = L2err(dtn,truescatt,dtn,dtnscatt,[dtn.theta, 2*pi]);
if verbose | abserr > 1e-12
    fprintf('L2 error for Rs = [1,2], V0s = [0,50]: %4.2e, %4.2e\n', abserr, relerr);
end

%% Test 6: two regions which make me sad, Rs = [2.5,6.5], V0s = [0,50]

% pick axisymmetric potential, piecewise constant on concentric annuli
Rs = [2.5,6.5]; V0s = [0,50]; 
N = length(V0s); coords = 'polar';
Vs = cell(1,N); for j = 1:N, Vs{j} = @(r,t) 0*r + V0s(j); end

% choose discretization parameters, make DtNBC object, solve
Nt = 40; Nrs = 0*V0s + 80;
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtnscatt = dtn.solve(k);

% figure out how many fourier coeffs we need
maxn = 40;
truescatt = compute_truescatt(k,V0s,dtn,maxn);
truescatt2 = compute_truescatt(k,V0s,dtn,2*maxn);
[abserr,relerr] = L2err(dtn,truescatt,dtn,truescatt2,[dtn.theta, 2*pi]);

part = @imag;
difference = truescatt-dtnscatt;
diffratio = difference./truescatt;

% collect maximum absolute and relative errors
[abserr,relerr] = L2err(dtn,truescatt,dtn,dtnscatt,[dtn.theta, 2*pi]);
if verbose | abserr > 1e-2
    fprintf('L2 error for Rs = [1,2], V0s = [0,50]: %4.2e, %4.2e (bad because inverting cond 1e-140 matrices?)\n', abserr, relerr);
    figure, dtn.plotValuesVec(log10(abs(part(difference))), 'final test absdiff',@real);
    figure, dtn.plotValuesVec(log10(abs(part(diffratio))), 'final test reldiff',@real);
end
% Since sqrt(k^2 - 50) ~ 6i, the argument to the bessel functions where 
% sqrt(k^2 - 50) and the outer radius are involved
% will be about 36i, where the bessel functions are large in magnitude.
% Then the matrix A will be extremely ill-conditioned.

warning on
disp('Warnings back on');

function test_chebBasis(verbose)

close all
addpath('..')

if nargin == 0, verbose = 0; end

%% f(x) = x^2
% cheb expansion is 0.5*T0(x) + 0.5*T2(x)
% fourier expansion is 0.25*exp(-2it) + 0.5 + 0.25*exp(2it), x = cos(t)
% derivative f'(x) = 2*x
% cheb expansion 2*T1(x)

f = @(x) x.^2;
c  = [0.5; 0; 0.5]; % nonzero part
dc = [0;2];

N = 5; % N + 1 points

[D,x] = cheb(N); % Chebyshev mesh
v = f(x); % values on full mesh
c_full = 0*v; c_full(1:length(c)) = c;
dc_full = 0*v; dc_full(1:length(dc)) = dc;
W = ones(length(x)); % map coeffs to values
W(:,2) = x;
for j = 3:length(x)
    W(:,j) = 2*x.*W(:,j-1) - W(:,j-2);
end

x2 = x((N+1)/2 + 1:end); % half mesh on [0,1]
v2 = f(x2); % values on half mesh
c_half = 0*x2; c_half([1,2]) = 0.5; % only even poly coeffs
 
%%
% Cheb coeffs from values on whole mesh

[abserr,relerr] = get_err(c_full, scattResComp2d.fullChebVals2chebCoeffs(v));
if verbose || abserr > 1e-15
    fprintf('--Error in cheb coeffs of f(x) = x^2 from full mesh: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(v, W*c_full);
if verbose || abserr > 1e-16
    fprintf('--Error in same comp with a matrix mult: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dc_full, W\(D*W*c_full));
if verbose || abserr > 1e-15
    fprintf('--Error in map from cheb coeffs to cheb coeffs of deriv: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% Cheb coeffs from values on half mesh plus symmetry condition

[abserr,relerr] = get_err(c_half, scattResComp2d.halfChebVals2chebCoeffs(v2,'even'));
if verbose || abserr > 1e-14
    fprintf('--Error in cheb coeffs of f(x) = x^2 from half mesh: %4.2e, %4.2e\n', abserr, relerr);
end

%% f(x) = x
% cheb expansion is f(x) = T1(x)
% fourier expansion is f(x) = 0.5*exp(-1i*t) + 0.5*exp(1i*t);

f = @(x) x;
c = [0; 1];

N = 5; % N + 1 points

[~,x] = cheb(N); % Chebyshev mesh
v = f(x); % values on full mesh
c_full = 0*v; c_full(1:length(c)) = c;

x2 = x((N+1)/2 + 1:end); % half mesh on [0,1]
v2 = f(x2); % values on half mesh
c_half = 0*x2; c_half(1) = 1; % only coeffs of odd polys

%%
% Cheb coeffs from values on whole mesh

[abserr,relerr] = get_err(c_full, scattResComp2d.fullChebVals2chebCoeffs(v));
if verbose || abserr > 1e-15
    fprintf('--Error in cheb coeffs of f(x) = x from full mesh: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% Cheb coeffs from values on half mesh plus symmetry condition

[abserr,relerr] = get_err(c_half, scattResComp2d.halfChebVals2chebCoeffs(v2,'odd'));
if verbose || abserr > 1e-14
    fprintf('--Error in cheb coeffs of f(x) = x from half mesh: %4.2e, %4.2e\n', abserr, relerr);
end

%% Both at once

f1 = @(x) x;    c1 = [0  ; 1; 0  ];
f2 = @(x) x.^2; c2 = [0.5; 0; 0.5];

N = 5;

[~,x] = cheb(N); % so six points
v = [f1(x), f2(x)];
c_full = 0*v; c_full(1:length(c1),:) = [c1, c2];

x2 = x(4:end);
v2 = v(4:end,:);
c_half = 0*v2; c_half(1:length(c1),:) = [c1, c2]; 

%%
% Values on whole mesh not using symmetry

[abserr,relerr] = get_err(c_full, scattResComp2d.fullChebVals2chebCoeffs(v));
if verbose || abserr > 1e-15
    fprintf('--Error in cheb coeffs of values matrix from full mesh: %4.2e, %4.2e\n', abserr, relerr);
end

%% f(x) = x^3 (antisymmetric) on piecewise radial meshes

f = @(x) x.^3;

%%
% first mesh

Nt = 6; Nrs = [5,11]; 
Rs = [1,1.75];
s1 = scattResComp2d(Nt,Nrs,Rs);
c1 = s1.radialVals2chebCoeffs(f(s1.r), 'odd');

%%
% second (bigger) mesh

Nrs = [7,13];
s2 = scattResComp2d(Nt,Nrs,Rs);
c2 = s2.radialVals2chebCoeffs(f(s2.r), 'odd');

%%
% compare

[abserr,relerr] = get_err(c1,c2([1:5, 8:end-2]));
if verbose || abserr > 1e-13
    fprintf('--Error in cheb coeff consistency of x^3 coeffs from two radial meshes: %4.2e, %4.2e\n', abserr, relerr);
end

%% Matrix with two columns

f = @(x) x.^3; g = @(x) x.^5; 

Nt = 6; Nrs = [7,11]; Rs = [0.5,1.75];
s1 = scattResComp2d(Nt,Nrs,Rs);
v1 = [f(s1.r), g(s1.r)];
c1 = s1.radialVals2chebCoeffs(v1, 'odd');

Nt = 6; Nrs = [8,13]; Rs = [0.5,1.75];
s2 = scattResComp2d(Nt,Nrs,Rs);
v2 = [f(s2.r), g(s2.r)];
c2 = s2.radialVals2chebCoeffs(v2, 'odd');

sub_idx = [1:7, 9:19];
[abserr,relerr] = get_err(c1, c2(sub_idx,:));
if verbose || abserr > 1e-13
    fprintf('--Error in cheb coeffs for x^3 and x^5 from different radial meshes: %4.2e, %4.2e\n', abserr, relerr);
end

%% Symmetric functions

f = @(x) 0*x+1; g = @(x) x.^2; h = @(x) x.^4;

Nt = 6; Nrs = [7,11]; Rs = [0.5,1.75];
s1 = scattResComp2d(Nt,Nrs,Rs);
v1 = [f(s1.r), g(s1.r), h(s1.r)];
c1 = s1.radialVals2chebCoeffs(v1, 'even');

Nt = 6; Nrs = [8,13]; Rs = [0.5,1.75];
s2 = scattResComp2d(Nt,Nrs,Rs);
v2 = [f(s2.r), g(s2.r), h(s2.r)];
c2 = s2.radialVals2chebCoeffs(v2, 'even');

sub_idx = [1:7, 9:19];
[abserr,relerr] = get_err(c1, c2(sub_idx,:));
if verbose || abserr > 1e-13
    fprintf('--Error in cheb coeffs for 1, x^2, x^4 from different radial meshes: %4.2e, %4.2e\n', abserr, relerr);
end

%% chebCoeffs2fullChebVals and fullChebVals2chebCoeffs

%% 
% even number of mesh points

v = rand(10,1);
c = scattResComp2d.fullChebVals2chebCoeffs(v);
w = scattResComp2d.chebCoeffs2fullChebVals(c);

[abserr,relerr] = get_err(v,w);
if verbose || abserr > 1e-14
    fprintf('--Error in consistency for fullChebVals2chebCoeffs and inverse (even mesh): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% odd number of mesh points

v = rand(11,1);
c = scattResComp2d.fullChebVals2chebCoeffs(v);
w = scattResComp2d.chebCoeffs2fullChebVals(c);

[abserr,relerr] = get_err(v,w);
if verbose || abserr > 1e-14
    fprintf('--Error in consistency for fullChebVals2chebCoeffs and inverse (odd mesh): %4.2e, %4.2e\n', abserr, relerr);
end

%% chebCoeffs2halfChebVals and halfChebVals2chebCoeffs

%% 
% even number of mesh points

N = 4; % size of half mesh
[~,x] = cheb(2*N-1); % 2N points on full mesh
x2 = x(N+1:end);
f = @(x) 5 + x.^2; % an even function
v = f(x2); % its values

c = scattResComp2d.halfChebVals2chebCoeffs(v,'even');
w = scattResComp2d.chebCoeffs2halfChebVals(c,'even');

[abserr,relerr] = get_err(v,w);
if verbose || abserr > 1e-14
    fprintf('--Error in consistency for halfChebVals2chebCoeffs and inverse (even mesh): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% odd number of mesh points

N = 5; % size of half mesh
[~,x] = cheb(2*N-1); % 2N points on full mesh
x2 = x(N+1:end);
f = @(x) 5 + x.^2; % an even function
v = f(x2); % its values

c = scattResComp2d.halfChebVals2chebCoeffs(v,'even');
w = scattResComp2d.chebCoeffs2halfChebVals(c,'even');

[abserr,relerr] = get_err(v,w);
if verbose || abserr > 1e-13
    fprintf('--Error in consistency for halfChebVals2chebCoeffs and inverse (odd mesh): %4.2e, %4.2e\n', abserr, relerr);
end

%% turn a fourierVec into its chebVec, axisymm function

f_even = @(x) cos(x) + 5*cos(x).^2;
f_odd = @(x) sin(x);

Nt = 6; Nrs = [20,20]; Rs = [0.5,1.75];
s1 = scattResComp2d(Nt,Nrs,Rs);
v1rect = zeros(s1.Nr,s1.Nt);
v1rect(:,2:2:end) = repmat(f_odd(s1.r), 1, s1.Nt/2);
v1rect(:,1:2:end) = repmat(f_even(s1.r), 1, s1.Nt/2);;
v1 = v1rect(:);
c1 = s1.chebVecFromFourierVec(v1);

Nt = 6; Nrs = [25,25]; Rs = [0.5,1.75];
s2 = scattResComp2d(Nt,Nrs,Rs);
v2rect = zeros(s2.Nr,s2.Nt);
v2rect(:,2:2:end) = repmat(f_odd(s2.r), 1, s2.Nt/2);
v2rect(:,1:2:end) = repmat(f_even(s2.r), 1, s2.Nt/2);
v2 = v2rect(:);
c2 = s2.chebVecFromFourierVec(v2);

sub_idx = [1:20, 26:45];
c2rect = reshape(c2,s2.Nr,s2.Nt);
c2rect = c2rect(sub_idx,:);

[abserr,relerr] = get_err(c1, c2rect(:));
if verbose || abserr > 1e-12
    fprintf('--Error for chebVecFromFourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

%% Consistency between chebVecFromFourierVec and fourierVecFromChebVec

[abserr,relerr] = get_err(v1, s1.fourierVecFromChebVec(c1));
if verbose || abserr > 1e-12
    fprintf('--Error in consistency between chebVecFromFourierVec and its inverse: %4.2e, %4.2e\n', abserr, relerr);
end

%% Consistency in first derivative

f = @(r,t) r.*exp(-1i*t) + r.^2 - r.*exp(1i*t) + 5*exp(2i*t);
df = @(r,t) exp(-1i*t) + 2*r - exp(1i*t);
coords = 'polar';

Nt = 4; Nrs = [4,5]; Rs = [0.8, 1.76];
s = scattResComp2d(Nt,Nrs,Rs);

% values
valuesVec  = s.valuesVecFromFun(f,coords);
fourierVec = s.fourierVecFromValuesVec(valuesVec);
chebVec    = s.chebVecFromFourierVec(fourierVec);

% check that the fourierVec is right
[abserr,relerr] = get_err(fourierVec, [s.r; s.r.^2; -s.r; 0*s.r+5]);
if verbose || abserr > 1e-14
    fprintf('--Error in recovering fourierVec from valuesVec: %4.2e, %4.2e\n', abserr, relerr);
end

% differentiate values with respect to r
DvaluesVec = s.Dr*valuesVec;
[abserr,relerr] = get_err(DvaluesVec, s.valuesVecFromFun(df,coords));
if verbose || abserr > 1e-11
    fprintf('--Error in Dr mapping values to derivatives: %4.2e, %4.2e\n', abserr, relerr);
end

% differentiate fourier with respect to r
DfourierRect = reshape(fourierVec,s.Nr,s.Nt);
DfourierRect(:,1:2:end) = s.Drs_odd *DfourierRect(:,1:2:end);
DfourierRect(:,2:2:end) = s.Drs_even*DfourierRect(:,2:2:end);
DfourierVec = DfourierRect(:);
[abserr,relerr] = get_err(DfourierVec, s.fourierVecFromValuesVec(DvaluesVec));
if verbose || abserr > 1e-12
    fprintf('--Error in mapping fourier values to derivatives: %4.2e, %4.2e\n', abserr, relerr);
end

% check derivative of fourierVec
DfourierVec_true = [0*s.r + 1; 2*s.r; 0*s.r - 1; 0*s.r];
[abserr,relerr] = get_err(DfourierVec,DfourierVec_true);
if verbose || abserr > 1e-13
    fprintf('--Error in DfourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

blehRect = reshape(chebVec,s.Nr,s.Nt);
% blehRect(:,2:2:end) = s.Winv_even*blehRect(:,2:2:end);
% blehRect(:,1:2:end) = s.Winv_odd *blehRect(:,1:2:end);
% [abserr,relerr] = get_err(fourierVec,blehRect(:))
blehRect(:,2:2:end) = s.Drs_even*s.Winv_even*blehRect(:,2:2:end);
blehRect(:,1:2:end) = s.Drs_odd *s.Winv_odd *blehRect(:,1:2:end);
[abserr,relerr] = get_err(DfourierVec,blehRect(:));
if verbose || abserr > 1e-13
    fprintf('--Error in how to use Drs_* and Winv_*: %4.2e, %4.2e\n', abserr, relerr);
end

DchebRect = reshape(chebVec,s.Nr,s.Nt);
DchebRect(:,1:2:end) = s.Drs_odd *s.Winv_odd *DchebRect(:,1:2:end);
DchebRect(:,2:2:end) = s.Drs_even*s.Winv_even*DchebRect(:,2:2:end);
DchebVec = s.chebVecFromFourierVec(DchebRect(:),1);
[abserr,relerr] = get_err(DchebVec, s.chebVecFromFourierVec(DfourierVec,1));
if verbose || abserr > 1e-13
    fprintf('--Error between DchebVec and converted DfourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

DchebRect = reshape(chebVec,s.Nr,s.Nt);
DchebRect(:,1:2:end) = s.Drs_cc_odd *DchebRect(:,1:2:end);
DchebRect(:,2:2:end) = s.Drs_cc_even*DchebRect(:,2:2:end);
DchebVec = DchebRect(:);
[abserr,relerr] = get_err(DchebVec,s.chebVecFromFourierVec(DfourierVec,1));
if verbose || abserr > 1e-13
    fprintf('--Error in using Drs_cc_*: %4.2e, %4.2e\n', abserr, relerr);
end

end

function [abserr,relerr] = get_err(x_true,x_approx)
    abserr = norm(full(x_true-x_approx));
    relerr = abserr/norm(full(x_true));
end

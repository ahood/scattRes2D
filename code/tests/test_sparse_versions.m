function test_sparse_versions(verbose)
% Want to be able to apply (A - E0*B) to a vector quickly. Why?
%   1) Eigenvalue computations. If E = E0 + 1/e is an eigenvalue, then
%   (A - E0*B)\B*x = e*x for x an eigenvector. Can use iterative method to
%   compute (A - E0*B)\y if we can do (A - E0*B)*y fast.
%   2) Scattering computations. Trying to solve (A - E*B)*x = b where b is
%   minus potential times incident wave stuff. x = (A - E*B)\b through
%   iterative method as above.

close all
addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end
disp('Temporarily turning warnings off');
warning off

%% Check U*x and apply_U(x) do the same thing

Nt = 10; Nrs = 20; Rs = rand;
scatt  = scattResComp2d(Nt,Nrs,Rs);
scatts = scattResComp2d_sparse(Nt,Nrs,Rs);

x = rand(scatt.Nt,1);
[abserr,relerr] = get_err(scatt.U*x,scatts.apply_U(x));
if verbose || abserr > 1e-12
    fprintf('Error between U and apply_U: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(scatt.Uinv*x, scatts.apply_Uinv(x));
if verbose || abserr > 1e-12
    fprintf('Error between Uinv and apply_Uinv: %4.2e, %4.2e\n', abserr, relerr);
end

%% Check same RHS is created
Nt = 10; Nrs = 20; Rs = rand;
Vs = {@(r,t) cos(r).*sin(t)};
coords = 'polar';
k = rand;
incfun = @(r,t) exp(1i*k*r.*cos(t));

% for dirBC and dirBC_sparse
dir = dirBC(Nt,Nrs,Vs,coords,Rs);
dirs = dirBC_sparse(Nt,Nrs,Vs,coords,Rs);

% Case 1: user passes wave number k instead of incfun
[abserr,relerr] = get_err(dir.RHSfromFun(k), dirs.RHSfromFun(k));
if verbose || abserr > 1e-10
    fprintf('Error for dirBC/dirBC_sparse RHSfromFun when k passed: %4.2e, %4.2e\n', abserr, relerr);
end

% Case 2: user passes incfun
[abserr,relerr] = get_err(dir.RHSfromFun(incfun), dirs.RHSfromFun(incfun));
if verbose || abserr > 1e-10
    fprintf('Error for dirBC/dirBC_sparse RHSfromFun when incfun passed: %4.2e, %4.2e\n', abserr, relerr);
end

% for pmlBC and pmlBC_sparse
Nr_out = 10;
Vs_pml = [Vs, {@(r,t) 0*r}];
Rout = Rs(end) + 1;
a = 0.4;
l = @(t) a*(t-Rout).^3;
dl = @(t) a*3*(t-Rout).^2;
d2l = @(t) a*3*2*(t-Rout);
pml = pmlBC(Nt,Nrs,Nr_out,Vs_pml,coords,Rs,Rout,l,dl,d2l);
pmls = pmlBC_sparse(Nt,Nrs,Nr_out,Vs_pml,coords,Rs,Rout,l,dl,d2l);

% Case 1: user passes wave number k instead of incfun
[abserr,relerr] = get_err(pml.RHSfromFun(k), pmls.RHSfromFun(k));
if verbose || abserr > 1e-10
    fprintf('Error for pmlBC/pmlBC_sparse RHSfromFun when k passed: %4.2e, %4.2e\n', abserr, relerr);
end

% Case 2: user passes incfun
[abserr,relerr] = get_err(pml.RHSfromFun(incfun), pmls.RHSfromFun(incfun));
if verbose || abserr > 1e-10
    fprintf('Error for pmlBC/pmlBC_sparse RHSfromFun when incfun passed: %4.2e, %4.2e\n', abserr, relerr);
end

% Case 3: user passes k and wants part of scattered wave on PML region too
[abserr,relerr] = get_err(pml.RHSfromFun_full(k), pmls.RHSfromFun_full(k));
if verbose || abserr > 1e-10
    fprintf('Error for pmlBC/pmlBC_sparse RHSfromFun_full when k passed: %4.2e, %4.2e\n', abserr, relerr);
end

% Case 4: user passes incfun and wants part of scattered wave on PML region
% too
[abserr,relerr] = get_err(pml.RHSfromFun_full(incfun), pmls.RHSfromFun_full(incfun));
if verbose || abserr > 1e-10
    fprintf('Error for pmlBC/pmlBC_sparse RHSfromFun_full when incfun passed: %4.2e, %4.2e\n', abserr, relerr);
end

% for DtNBC and DtNBC_sparse
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dtns = DtNBC_sparse(Nt,Nrs,Vs,coords,Rs);

% Case 1: user passes wave number k instead of incfun
[abserr,relerr] = get_err(dtn.RHSfromFun(k), dtns.RHSfromFun(k));
if verbose || abserr > 1e-10
    fprintf('Error for DtNBC/DtNBC_sparse RHSfromFun when k passed: %4.2e, %4.2e\n', abserr, relerr);
end

% Case 2: user passes incfun
[abserr,relerr] = get_err(dtn.RHSfromFun(incfun), dtns.RHSfromFun(incfun));
if verbose || abserr > 1e-10
    fprintf('Error for DtNBC/DtNBC_sparse RHSfromFun when incfun passed: %4.2e, %4.2e\n', abserr, relerr);
end

%% Check solve functions return same results

[abserr,relerr] = get_err(dir.solve(k),dirs.solve(k));
if verbose || abserr > 1e-10
    fprintf('Error between dir.solve(k) and dirs.solve(k): %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(pml.solve_full(k),pmls.solve_full(k));
if verbose || abserr > 1e-10
    fprintf('Error between pml.solve_full(k) and pmls.solve_full(k): %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dtn.solve(k),dtns.solve(k));
if verbose || abserr > 1e-10
    fprintf('Error between dtn.solve(k) and dtns.solve(k): %4.2e, %4.2e\n', abserr, relerr);
end

%% Check converting between representations behaves the same
% Testing
% 1. valuesVecFromFourierVec and fourierVecFromValuesVec
% 2. valuesVecFromFun and fourierVecFromFun
% 3. valuesVecFromFunCellArray and fourierVecFromFunCellArray
% 4. chebVecFromFourierVec and fourierVecFromChebVec

% set up instances
Nt = 4; Nrs = [20,25]; Rs = [1.1,2];
s = scattResComp2d(Nt,Nrs,Rs);
s_sparse = scattResComp2d(Nt,Nrs,Rs);

% pick a function for evaluation, defined first in term of its fourier
% coefficient functions, which are even functions for even index and odd
% for odd index
fm1 = @(r) r; 
f0  = @(r) r.^2 + 5;
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
[abserr,relerr] = get_err(s.valuesVecFromFourierVec(exactFourierVec),s_sparse.valuesVecFromFourierVec(exactFourierVec));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(s.fourierVecFromValuesVec(exactValuesVec),s_sparse.fourierVecFromValuesVec(exactValuesVec));
if verbose || abserr > 1e-16
    fprintf('Error for fourierVecFromValuesVec: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% 2. valuesVecFromFun and fourierVecFromFun
[abserr,relerr] = get_err(s.valuesVecFromFun(f,'polar'),s_sparse.valuesVecFromFun(f,'polar'));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFun (polar):   %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(s.valuesVecFromFun(frect,'rect'),s_sparse.valuesVecFromFun(frect,'rect'));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFun (rect):    %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(s.valuesVecFromFun(fcomp,'complex'),s_sparse.valuesVecFromFun(fcomp,'complex'));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFun (complex): %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(s.fourierVecFromFun(f,'polar'),s_sparse.fourierVecFromFun(f,'polar'));
if verbose || abserr > 1e-16
    fprintf('Error for fourierVecFromFun (polar):  %4.2e, %4.2e\n', abserr, relerr);
end

%%
% 3. valuesVecFromFunCellArray and fourierVecFromFunCellArray
[abserr,relerr] = get_err(s.valuesVecFromFunCellArray(fs,'polar'),s_sparse.valuesVecFromFunCellArray(fs,'polar'));
if verbose || abserr > 1e-16
    fprintf('Error for valuesVecFromFunCellArray (polar): %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(s.fourierVecFromFunCellArray(fs,'polar'),s_sparse.fourierVecFromFunCellArray(fs,'polar'));
if verbose || abserr > 1e-16
    fprintf('Error for fourierVecFromFunCellArray (polar): %4.2e, %4.2e\n', abserr, relerr);
end

%%
% 4. chebVecFromFourierVec and fourierVecFromChebVec

chebVec = s.chebVecFromFourierVec(exactFourierVec);

[abserr,relerr] = get_err(chebVec, s_sparse.chebVecFromFourierVec(exactFourierVec));
if verbose || abserr > 1e-16
    fprintf('Error for chebVecFromFourierVec: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(s.fourierVecFromChebVec(exactFourierVec), ...
                          s_sparse.fourierVecFromChebVec(exactFourierVec));
if verbose || abserr > 1e-16
    fprintf('Error for fourierVecFromChebVec: %4.2e, %4.2e\n', abserr, relerr);
end

% [abserr,relerr] = get_err(s.Winv, s_sparse.Winv);
% if verbose || abserr > 1e-16
%     fprintf('Error in Winv consistency: %4.2e, %4.2e\n', abserr, relerr);
% end

%%
% choose piecewise constant, axisymmetric 2D potential
Rin    = 1.5; 
Rout   = 2.5; 
V0_in  = 0; 
V0_out = 38;
coords = 'polar';
Vs = {@(r,t) 0*r + V0_in; @(r,t) 0*r + V0_out};
pot_params = struct('Vs',{Vs},'Rs',[Rin,Rout],'Vvals',[V0_in,V0_out]);

% Choosing discretization of potential support
% Nt = 40;
% Nrs = [20,20]; 
Nt = 6;
Nrs = [4, 4];

% DtN problem
dtn = DtNBC(Nt,Nrs,pot_params.Vs,coords,pot_params.Rs);

% sparse DtN problem
dtns = DtNBC_sparse(Nt,Nrs,pot_params.Vs,coords,pot_params.Rs);

% random vector to multiply, random frequency parameter
x0 = rand(dtn.Nt*dtn.Nr,1);
k = 2*rand + 1;
% k = 2;
z = k^2; % need this for setting up rat approx

% apply T(k) two ways and compare results
if verbose
    disp('----------------------------------------------------------------')
    disp('Testing DtN multiplication...dtn.T(k)*x0 vs. dtns.apply_T(x0,k)')
    disp('----------------------------------------------------------------')
end
x_dtn  = dtn.T(k)*x0;
x_dtns = dtns.apply_T(x0,k);

abserr = norm(x_dtn - x_dtns); relerr = abserr/norm(x_dtn);
if verbose || abserr > 1e-12
    fprintf('k is %f\n', k);
    disp('dtn.T(k)*x0 vs. dtns.apply_T(x0,k)');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

if verbose
    disp('----------------------------------------------------------------')
    disp('Testing DtN inversion...dtn.T(k)\x0 vs. gmres with dtns.apply_T')
    disp('----------------------------------------------------------------')
end
x_dtn  = dtn.T(k)\x0;
x_dtns = dtns.solve(k,x0,[],verbose); % default preconditioner

abserr = norm(x_dtn - x_dtns); relerr = abserr/norm(x_dtn);
if verbose || abserr > 1e-5
    fprintf('k is %f\n', k);   
    disp('dtn.T(k)\x0 vs. gmres with dtns.apply_T');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

% do same for pml stuff

% choose discretization parameters and PML parameters
RPML = 8; 
Nr_out = 5;
a = 0.4;
l = @(t) a*(t-Rout).^3;
dl = @(t) a*3*(t-Rout).^2;
d2l = @(t) a*3*2*(t-Rout);
pmlVs = Vs; pmlVs{end+1} = @(r,t) 0*r; % isn't this wrong?
pml_params = struct('pmlVs',{pmlVs},'l',l,'dl',dl,'d2l',d2l,...
    'Nr_out',Nr_out,'RPML',RPML);
if verbose
    disp('----------------------------------------------------------------')
    disp('Testing PML multiplication...pml.Tfull(k)*x0 vs. pmls.apply_Tfull(x0,k)')
    disp('----------------------------------------------------------------')
end
pml  = pmlBC       (Nt,dtn.Nrs,Nr_out,Vs,coords,dtn.Rs,RPML,l,dl,d2l);
pmls = pmlBC_sparse(Nt,dtn.Nrs,Nr_out,Vs,coords,dtn.Rs,RPML,l,dl,d2l);
x0 = rand(pml.Nr*pml.Nt,1);
x_pml = pml.Tfull(k)*x0;
x_pmls = pmls.apply_Tfull(x0,k);

abserr = norm(x_pml - x_pmls); relerr = abserr/norm(x_pml);
if verbose || abserr > 1e-8
    fprintf('k is %f\n', k);
    disp('pml.Tfull(k)*x0 vs. pmls.apply_Tfull(x0,k)');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

if verbose
    disp('----------------------------------------------------------------')
    disp('Testing PML inversion...pml.Tfull(k)\x0 vs. gmres with pmls.apply_T')
    disp('----------------------------------------------------------------')
end
x_pml = pml.Tfull(k)\x0;
x_pmls = pmls.solve_full(k,x0,[],verbose); % default preconditioner

abserr = norm(x_pml - x_pmls); relerr = abserr/norm(x_pml);
if verbose || abserr > 1e-4
    fprintf('k is %f\n', k);
    disp('pml.Tfull(k)\x0 vs. gmres with pmls.apply_T');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

% do same for dir stuff
dir  = dirBC       (Nt,Nrs,Vs,coords,dtn.Rs);
dirs = dirBC_sparse(Nt,Nrs,Vs,coords,dtn.Rs);
x0 = rand(dir.Nr*dir.Nt,1);

if verbose
    disp('----------------------------------------------------------------')
    disp('Testing dir multiplication...dir.T(k)*x0 vs. dirs.apply_T(x0,k)')
    disp('----------------------------------------------------------------')
end
x_dir = dir.T(k)*x0;
x_dirs = dirs.apply_T(x0,k);

abserr = norm(x_dir - x_dirs); relerr = abserr/norm(x_dir);
if verbose || abserr > 1e-12
    fprintf('k is %f\n', k);    
    disp('dir.T(k)*x0 vs. dirs.apply_T(x0,k)');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

if verbose
    disp('----------------------------------------------------------------')
    disp('Testing dir inversion...dir.T(k)\x0 vs. gmres with dirs.apply_T')
    disp('----------------------------------------------------------------')
end
x_dir  = dir.T(k)\x0;
x_dirs = dirs.solve(k,x0,[],verbose); % default preconditioner

abserr = norm(x_dir - x_dirs); relerr = abserr/norm(x_dir);
if verbose || abserr > 1e-5
    fprintf('k is %f\n', k);
    disp('dir.T(k)\x0 vs. gmres with dirs.apply_T');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

% do same for rat approx to DtN map
ell = ellipse(z,0,1.5,1.5,10,10);
br = 1;

N = 40; % small number of points for rational approxes
rat  = ratApproxDtNBC       (dtn ,ell,br,N);
rats = ratApproxDtNBC_sparse(dtns,ell,br,N);

if verbose
    disp('----------------------------------------------------------------')
    disp('Testing rat multiplication...rat.T(k)*x0 vs. dtn.T(k)*x0')
    disp('----------------------------------------------------------------')
end
x_rat = rat.T(k)*x0;
x_dtn = dtn.T(k)*x0;
abserr = norm(x_rat - x_dtn); relerr = abserr/norm(x_dtn);
if verbose || abserr > 1e-10
    fprintf('k is %f\n', k);
    disp('rat.T(k)*x0 vs. dtn.T(k)*x0');
    fprintf('  Error between rat and dtn: %4.2e, %4.2e\n', abserr, relerr);
end

if verbose
    disp('----------------------------------------------------------------')
    disp('Testing rat multiplication...rat.T(k)*x0 vs. rats.apply_T(x0,k)')
    disp('----------------------------------------------------------------')
end
x_rats = rats.apply_T(x0,k);
abserr = norm(x_rat - x_rats); relerr = abserr/norm(x_rat);
if verbose || abserr > 1e-13
    fprintf('k is %f\n', k);
    disp('rat.T(k)*x0 vs. rats.apply_T(x0,k)');
    fprintf('  Error between rat and rats: %4.2e, %4.2e\n', abserr, relerr);
end

% rat.T(k)\x0 vs rat.solve(k,x0)
if verbose
    disp('----------------------------------------------------------------')
    disp('Testing rat scattering...rat.T(k)\x0 vs. rats.solve(k,x0)       ')
    disp('----------------------------------------------------------------')
end

x_rat = rat.T(k)\x0;
x_rats = rats.solve(k,x0);
abserr = norm(x_rat - x_rats); relerr = abserr/norm(x_rat);
if verbose || abserr > 1e-13
    fprintf('k is %f\n', k);
    disp('rat.T(k)*x0 vs. rats.apply_T(x0,k)');
    fprintf('  Error between rat and rats: %4.2e, %4.2e\n', abserr, relerr);
end

% do same for matrix whose schur comp gives rat approx to DtN map
if verbose
    disp('----------------------------------------------------------------------------')
    disp('Testing rat full multiplication...rat.Tfull(k)*x0 vs. rats.apply_Tfull(x0,k)')
    disp('----------------------------------------------------------------------------')
end
x0 = rand(length(rat.A),1);

x_rat = rat.Tfull(k)*x0;
x_rats = rats.apply_Tfull(x0,k);

abserr = norm(x_rat - x_rats); relerr = abserr/norm(x_rat);
if verbose || abserr > 1e-13
    fprintf('k is %f\n', k);
    disp('rat.Tfull(k)*x0 vs. rats.apply_Tfull(x0,k)');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

if verbose
    disp('-----------------------------------------------------------------------')
    disp('Testing rat inversion...rat.Tfull(k)\x0 vs. gmres with rats.apply_Tfull')
    disp('-----------------------------------------------------------------------')
end
x_rat = rat.Tfull(k)\x0;
Afun = @(x) rats.apply_Tfull(x,k);
if verbose
    x_rats = gmres(Afun,x0,[],[],1000); 
else
    [x_rats,flag] = gmres(Afun,x0,[],[],1000); 
end

abserr = norm(x_rat - x_rats); relerr = abserr/norm(x_rat);
if verbose || abserr > 1e-5
    fprintf('k is %f\n', k);
    disp('rat.Tfull(k)\x0 vs. gmres with rats.apply_Tfull');
    fprintf('  Error: %4.2e, %4.2e\n', abserr, relerr);
end

if verbose
    disp('--------------------------------------------------')
    disp('Testing Dirichlet eigenvalue computations         ')
    disp('--------------------------------------------------')
end

E0 = ell.c;
n = 3;

% compute all eigs and get the n closest to E0
dir.eig_comp();
dir_eigs = dir.get_closest(dir.evals,n,E0);
dir_eigs = sort(dir_eigs);

% same but iteratively -- have to sort
[~,dirs_eigs] = dirs.eig_comp(n,E0,0); % 10 eigs closest to 2, suppress GMRES output
dirs_eigs = sort(dirs_eigs);

abserr = norm(dir_eigs - dirs_eigs);
relerr = abserr/norm(dir_eigs);
if verbose || abserr > 1e-4
    [dir_eigs.'; dirs_eigs.']
    fprintf('  Error between dirBC and dirBC_sparse eigs: %4.2e, %4.2e\n', abserr,relerr);
end

if verbose
    disp('--------------------------------------------')
    disp('Testing PML eigenvalue computations         ')
    disp('--------------------------------------------')
end

% compute all eigs and get the n closest to E0
pml.eig_comp();
pml_eigs = pml.get_closest(pml.evals,n,E0);
pml_eigs = sort(pml_eigs);

% same but iteratively -- have to sort
[~,pmls_eigs] = pmls.eig_comp(n,E0,0); % 10 eigs closest to 2, suppress GMRES output
pmls_eigs = sort(pmls_eigs);

abserr = norm(pml_eigs - pmls_eigs);
relerr = abserr/norm(pml_eigs);
if verbose || abserr > 1e-4
    [pml_eigs.'; pmls_eigs.']
    fprintf('  Error between pmlBC and pmlBC_sparse eigs: %4.2e, %4.2e\n', abserr,relerr);
end

% compute the n resonance energies (k^2) closest to ell.c
same_number_converged = false;
while ~same_number_converged
    % by feeding A and B to eig, where rat.Tfull(k) = A - k^2*B    
    tic; 
    if verbose
        fprintf('--Computing %d eigs closest to ellipse center using formed matrix...',n);
    end
    opts.isreal = false;
    [~,rat_res] = rat.resonances(n);

    if verbose
        fprintf('  --Took %.1f s to get evals\n', toc);   
    end

    tic; 
    if verbose
        fprintf('--Computing %d eigs closest to ellipse center using gmres...\n',n);
    end

    [~,rats_res] = rats.resonances(n);

    if verbose
        fprintf('  --Took %.1f s to get evals\n', toc);   
    end

    % compare
    rat_res  = sort(rat_res ( ~isnan(rat_res) ));
    rats_res = sort(rats_res( ~isnan(rats_res) ));

    if length(rat_res) == length(rats_res)
        same_number_converged = true;
    else
        if verbose
            disp('    **Different numbers of eigenvalues converged. Trying again.');
        end
    end
end

abserr = norm(rat_res - rats_res); relerr = abserr/norm(rat_res);
if verbose || abserr > 1e-4
    [rat_res.'; rats_res.']
    fprintf('  --Error between the two sets of eigs: %4.2e, %4.2e\n', abserr, relerr);
end

warning on
disp('Warnings back on');

end

function x = Tfull_k0_inv_fun(Tfull_k0,x,itmax,M_k0) 
    [x,flag] = gmres(Tfull_k0, x,[],1e-6,itmax,M_k0);
end

function [abserr,relerr] = get_err(x_true,x_approx)
    abserr = norm(full(x_true-x_approx));
    relerr = abserr/norm(full(x_true));
end
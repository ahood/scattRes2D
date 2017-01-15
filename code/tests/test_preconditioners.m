function test_preconditioners(verbose)

close all
addpath(fileparts(pwd));
if nargin == 0, verbose = 0; end

%% Piecewise constant axisymmetric potential

if verbose
    disp('Testing piecewise constant axisymmetric potential...')
end

Rin    = 1.5; 
Rout   = 2.5; 
V0_in  = 0; 
V0_out = 38;
coords = 'polar';
Vs = {@(r,t) 0*r + V0_in; @(r,t) 0*r + V0_out};
pot_params = struct('Vs',{Vs},'Rs',[Rin,Rout],'Vvals',[V0_in,V0_out]);

% Choosing discretization of potential support
Nt = 10; Nrs = [5,5];

%%
% DtN problems

dtn = DtNBC(Nt,Nrs,pot_params.Vs,coords,pot_params.Rs);
dtns = DtNBC_sparse(Nt,Nrs,pot_params.Vs,coords,pot_params.Rs);

%%
% random frequency parameter
k = 2*rand + 1;
z = k^2;
if verbose
    fprintf('Frequency parameter k = %.6f\n', k);
end

%%
% rational problems

% the ellipse where we put poles for rat approx to DtN map
ell = ellipse(z,0,1,1,10,10); % circle centered at z
br = 1;

rat  = ratApproxDtNBC       (dtn ,ell,br);
rats = ratApproxDtNBC_sparse(dtns,ell,br);

%%
% random vectors which will solve A*x = b (b will be computed in a minute)
x0 = rand(dtns.Nr*dtns.Nt,1); 
x0_full = [x0; rand(length(rats.A22diag),1)];

%% 
% apply dtn.T(k) two ways and compare results

b_dtn  = dtn.T(k)*x0;
b_dtns = dtns.apply_T(x0,k);

abserr = norm(b_dtn - b_dtns);
relerr = abserr/norm(b_dtn);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)*x0 and dtns.apply_T(x0,k): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve dtn.T(k)*x = b with block diagonal preconditioner

x_dtns = dtns.solve(k,b_dtn,@(x) dtns.apply_pc_blkdiag(x,k),verbose);

abserr = norm(x_dtns - x0);
relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)\\b and iterative sol with blkdiag precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve dtn.T(k)*x = b with fourier preconditioner, av potential, and the part of boundary conditions that looks block diagonal in Fourier space

% compare action to that of T(k)
Tx = dtns.apply_T(x0,k);
Mx = dtns.apply_pc_fourier_Vav_blkdiagICBC(x0,k,'mult');

abserr = norm(Mx-Tx); relerr = abserr/norm(Tx);
if verbose || abserr > 1e-6
    fprintf('--Error between mult action of dtn.T and fourier precond (surprisingly small!): %4.2e, %4.2e\n', abserr, relerr);
end

% compare action of M(k)^{-1} to that of T(k)^{-1}
Minvb = dtns.apply_pc_fourier_Vav_blkdiagICBC(b_dtn,k,'inv');

abserr = norm(Minvb - x0); relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)\\b and (fourier precond)\\b (hella awesome!): %4.2e, %4.2e\n', abserr, relerr);
end

% solve with preconditioner
x_dtns = dtns.solve(k,b_dtn,@(x) dtns.apply_pc_fourier_Vav_blkdiagICBC(x,k,'inv'),verbose);

abserr = norm(x_dtns - x0);
relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)\\b and iterative sol with fourier precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% apply rat.T(k) two ways and compare results

b_rat  = rat.T(k)*x0;
b_rats = rats.apply_T(x0,k);

abserr = norm(b_rat - b_rats);
relerr = abserr/norm(b_rat);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.T(k)*x0 and rats.apply_T(x0,k): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.T(k)*x = x0 with DtN block diagonal preconditioner

x_rats = rats.solve(k,b_rat,@(x) dtns.apply_pc_blkdiag(x,k), verbose);

abserr = norm(x0 - x_rats);
relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.T(k)\\b and iterative sol with DtN blkdiag precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.T(k)*x = x0 with DtNBC fourier preconditioner, av potential, and the part of boundary conditions that looks block diagonal in Fourier space

x_rats = rats.solve(k,b_rat,@(x) rats.apply_scpc_fourier_Vav_blkdiagICBC(x,k,'inv'), verbose);

abserr = norm(x_rats - x0); relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.T(k)\\b and iterative sol with Fourier precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% apply rat.Tfull(k) two ways and compare results

b_rat  = rat.Tfull(k)*x0_full;
b_rats = rats.apply_Tfull(x0_full,k);

abserr = norm(b_rat - b_rats);
relerr = abserr/norm(b_rat);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.Tfull(k)*x0 and rats.apply_Tfull(x0,k): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.Tfull(k)*x = b with block diagonal preconditioner

% no solve_full because I don't know how to make sense of the right hand
% side physically
itmax = length(b_rat);
if verbose
    x_rats = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_blkdiag(x,k));
else
    [x_rats,~] = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_blkdiag(x,k));
end

abserr = norm(x0_full - x_rats);
relerr = abserr/norm(x0_full);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.Tfull(k)\\b and iterative sol with blkdiag precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.T_full(k)*x = x0 with fourier preconditioner, av potential, and the part of boundary conditions that looks block diagonal in Fourier space

itmax = length(b_rat);
if verbose
    x_rats = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_fourier_Vav_blkdiagICBC(x,k,'inv'));
else
    [x_rats,~] = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_fourier_Vav_blkdiagICBC(x,k,'inv'));
end

abserr = norm(x_rats - x0_full); relerr = abserr/norm(x0_full);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.Tfull(k)\\b and M(k)\\b with fourier precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% Three bump potential

if verbose
    disp('Testing three bump potential from Lin''s paper...')
end

% potential from Lin's paper
Rcenters = 1.4;
sigma = 1/3;
hbar = 0.025;
laplaceCoeff = hbar^2/2;
c = Rcenters*exp(2i*pi*(1:3)/3);
G = @(x,y,cj) exp(-abs(x+1i*y-cj).^2/2/sigma^2);
V = @(x,y) (G(x,y,c(1)) + G(x,y,c(2)) + G(x,y,c(3)))/laplaceCoeff;
coords = 'rect';
cdist = 2.5;
R = Rcenters + cdist;
Nt = 10; Nr = 10; 

%%
% DtN problems

dtn = DtNBC(Nt,Nr,{V},coords,R);
dtns = DtNBC_sparse(Nt,Nr,{V},coords,R);

%%
% typical frequency parameter
k = sqrt(1600-40i);
z = k^2;
if verbose
    fprintf('Frequency parameter k = %.6f\n', k);
end

%%
% rational problems

% the ellipse where we put poles for rat approx to DtN map
ell = ellipse(z,0,1,1,10,10); % circle centered at z
br = 1;

rat  = ratApproxDtNBC       (dtn ,ell,br);
rats = ratApproxDtNBC_sparse(dtns,ell,br);

%%
% random vectors for testing multiplication and inversion actions
x0 = rand(dtns.Nr*dtns.Nt,1);
x0_full = [x0; rand(length(rats.A22diag),1)];

%% 
% apply dtn.T(k) two ways and compare results

b_dtn  = dtn.T(k)*x0;
b_dtns = dtns.apply_T(x0,k);

abserr = norm(b_dtn - b_dtns);
relerr = abserr/norm(b_dtn);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)*x0 and dtns.apply_T(x0,k): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve dtn.T(k)*x = b with block diagonal preconditioner

x_dtns = dtns.solve(k,b_dtn,@(x) dtns.apply_pc_blkdiag(x,k),verbose);

abserr = norm(x_dtns - x0);
relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)\\b and iterative sol with blkdiag precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve dtn.T(k)*x = b with fourier preconditioner, av potential, and the part of boundary conditions that looks block diagonal in Fourier space

x_dtns = dtns.solve(k,b_dtn,@(x) dtns.apply_pc_fourier_Vav_blkdiagICBC(x,k,'inv'),verbose);

abserr = norm(x_dtns - x0);
relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between dtn.T(k)\\b and iterative sol with fourier precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% apply rat.T(k) two ways and compare results

b_rat  = rat.T(k)*x0;
b_rats = rats.apply_T(x0,k);

abserr = norm(b_rat - b_rats);
relerr = abserr/norm(b_rat);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.T(k)*x0 and rats.apply_T(x0,k): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.T(k)*x = x0 with block diagonal preconditioner

x_rats = rats.solve(k,b_rat,@(x) dtns.apply_pc_blkdiag(x,k),verbose);

abserr = norm(x0 - x_rats);
relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.T(k)\\b and iterative sol with blkdiag precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.T(k)*x = x0 with DtNBC fourier preconditioner, av potential, and the part of boundary conditions that looks block diagonal in Fourier space

x_rats = rats.solve(k,b_rat,@(x) rats.apply_scpc_fourier_Vav_blkdiagICBC(x,k,'inv'),verbose);

abserr = norm(x_rats - x0); relerr = abserr/norm(x0);
if verbose || abserr > 1e-6
    fprintf('--Error between T(k)\\b and iterative sol with fourier precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% apply rat.Tfull(k) two ways and compare results

b_rat  = rat.Tfull(k)*x0_full;
b_rats = rats.apply_Tfull(x0_full,k);

abserr = norm(b_rat - b_rats);
relerr = abserr/norm(b_rat);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.Tfull(k)*x0 and rats.apply_Tfull(x0,k): %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.Tfull(k)*x = b with block diagonal preconditioner

itmax = length(b_rat);
if verbose
    x_rats = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_blkdiag(x,k));
else
    [x_rats,~] = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_blkdiag(x,k));
end

abserr = norm(x0_full - x_rats);
relerr = abserr/norm(x0_full);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.Tfull(k)\\b and iterative sol with blkdiag precond: %4.2e, %4.2e\n', abserr, relerr);
end

%% 
% Iteratively solve rat.T_full(k)*x = x0 with fourier preconditioner, av potential, and the part of boundary conditions that looks block diagonal in Fourier space

itmax = length(b_rat);
if verbose
    x_rats = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_fourier_Vav_blkdiagICBC(x,k,'inv'));
else
    [x_rats,~] = gmres(@(x) rats.apply_Tfull(x,k),b_rat,[],1e-12,itmax,@(x) rats.apply_pc_fourier_Vav_blkdiagICBC(x,k,'inv'));
end

abserr = norm(x_rats - x0_full); relerr = abserr/norm(x0_full);
if verbose || abserr > 1e-6
    fprintf('--Error between rat.Tfull(k)\\b and iterative sol with fourier precond: %4.2e, %4.2e\n', abserr, relerr);
end


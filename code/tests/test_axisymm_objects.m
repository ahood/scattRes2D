function test_axisymm_objects(verbose)

addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end
disp('Temporarily turning warnings off');
warning off

%% More general example
if verbose, disp('More general example...'); end

% discretization parameters
Nt = 50;
Nrs = 40;

% axisymmetric 2D potential
Rs = 3.9;
coords = 'polar';
Vs = {@(r,t) r.^2 - 3*r};

% PML parameters
Rout = Rs(end); RPML = 8; Nr_out = 40; a = 0.4;
l = @(t) a*(t-Rout).^3; dl = @(t) a*3*(t-Rout).^2; d2l = @(t) a*3*2*(t-Rout);
pmlVs = Vs; pmlVs{end+1} = @(r,t) 0*r;

% pick energy
k = 2;

%% Dir BCs

dir         = dirBC        (Nt,Nrs,Vs,coords,Rs);
dir_axisymm = dirBC_axisymm(Nt,Nrs,Vs,coords,Rs);

%%
% In Fourier basis

dir_scattFourierVec         = dir.solve_fc(k);
dir_scattFourierVec_axisymm = dir_axisymm.solve_fc(k);

[abserr,relerr] = get_err(dir_scattFourierVec, dir_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between dirBC and dirBC_axisymm sols in Fourier basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

dir_scattChebVec         = dir.solve_cc(k);
dir_scattChebVec_axisymm = dir_axisymm.solve_cc(k);

[abserr,relerr] = get_err(dir_scattChebVec, dir.chebVecFromFourierVec(dir_scattFourierVec));
if verbose || abserr > 1e-5
    fprintf('Error in consistency between sols in Fourier and Cheb bases with dir BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dir_scattChebVec, dir_scattChebVec_axisymm);
if verbose || abserr > 1e-12
    fprintf('Error between dirBC and dirBC_axisymm sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

%% PML with dir BCs

pml         = pmlBC        (Nt,Nrs,Nr_out,pmlVs,coords,Rs,RPML,l,dl,d2l);
pml_axisymm = pmlBC_axisymm(Nt,Nrs,Nr_out,pmlVs,coords,Rs,RPML,l,dl,d2l);

%%
% In Fourier basis

pmlfull_scattFourierVec         = pml.solve_full_fc(k);
pmlfull_scattFourierVec_axisymm = pml_axisymm.solve_full_fc(k);

[abserr,relerr] = get_err(pmlfull_scattFourierVec, pmlfull_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm full sols in Fourier basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

pmlfull_scattChebVec         = pml.solve_full_cc(k);
pmlfull_scattChebVec_axisymm = pml_axisymm.solve_full_cc(k);

[abserr,relerr] = get_err(pmlfull_scattChebVec, pmlfull_scattChebVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm full sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

%% PML BC from Schur complement

%%
% In Fourier basis

pml_scattFourierVec         = pml.solve_fc(k);
pml_scattFourierVec_axisymm = pml_axisymm.solve_fc(k);
pml_scattValuesVec = pml.solve(k);

[abserr,relerr] = get_err(pml_scattFourierVec, pml_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm sols in Fourier basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

pml_scattChebVec         = pml.solve_cc(k);
pml_scattChebVec_axisymm = pml_axisymm.solve_cc(k);

[abserr,relerr] = get_err(pml_scattChebVec, dir.chebVecFromFourierVec(pml_scattFourierVec));
if verbose || abserr > 1e-5
    fprintf('Error in consistency between sols in Fourier and Cheb bases with pml BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(pml_scattChebVec, pml_scattChebVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

%% DtN BCs

dtn         = DtNBC        (Nt,Nrs,Vs,coords,Rs);
dtn_axisymm = DtNBC_axisymm(Nt,Nrs,Vs,coords,Rs);

%%
% In Fourier basis

dtn_scattFourierVec         = dtn.solve_fc(k);
dtn_scattFourierVec_axisymm = dtn_axisymm.solve_fc(k);
dtn_scattValuesVec = dtn.solve(k);

[abserr,relerr] = get_err(dtn_scattFourierVec, dtn_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between DtNBC and DtNBC_axisymm sols: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = get_err(dtn_scattValuesVec, pml_scattValuesVec);
if verbose || abserr > 1e-6
    fprintf('Error between DtNBC and pmlBC sols in values basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

dtn_scattChebVec         = dtn.solve_cc(k);
dtn_scattChebVec_axisymm = dtn_axisymm.solve_cc(k);

[abserr,relerr] = get_err(dtn_scattChebVec, dtn.chebVecFromFourierVec(dtn_scattFourierVec));
if verbose || abserr > 1e-5
    fprintf('Error in consistency between sols in Fourier and Cheb bases with DtN BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dtn_scattChebVec, dtn_scattChebVec_axisymm);
if verbose || abserr > 1e-12
    fprintf('Error between DtNBC and DtNBC_axisymm sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dtn.C(k), pml.C(k));
if verbose || abserr > 1e-5
    fprintf('Error between dtn.C and pml.C at k = %f: %4.2e, %4.2e\n', k, abserr, relerr);
end

%% Piecewise constant example

if verbose, disp('piecewise constant example...'); end

% discretization parameters
Nt = 50;
Nrs = [20,20];

% piecewise constant, axisymmetric 2D potential
Rin = 1.4; Rout = 3.9; V0_in = 0; V0_out = 38;
Rs = [Rin,Rout];
coords = 'polar';
Vs = {@(r,t) 0*r + V0_in; @(r,t) 0*r + V0_out};

% PML parameters
RPML = 8; Nr_out = 40; a = 0.4;
l = @(t) a*(t-Rout).^3; dl = @(t) a*3*(t-Rout).^2; d2l = @(t) a*3*2*(t-Rout);
pmlVs = Vs; pmlVs{end+1} = @(r,t) 0*r;

% pick energy and make true solution
k = 2;
s = scattResComp2d(Nt,Nrs,Rs);
true_scattValuesVec = compute_truescatt(k,[V0_in,V0_out],s,Nt/2);

%% Dir BCs

dir         = dirBC        (Nt,Nrs,Vs,coords,Rs);
dir_axisymm = dirBC_axisymm(Nt,Nrs,Vs,coords,Rs);

%%
% In Fourier basis

dir_scattFourierVec         = dir.solve_fc(k);
dir_scattFourierVec_axisymm = dir_axisymm.solve_fc(k);

[abserr,relerr] = get_err(dir_scattFourierVec, dir_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between dirBC and dirBC_axisymm sols in Fourier basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

dir_scattChebVec         = dir.solve_cc(k);
dir_scattChebVec_axisymm = dir_axisymm.solve_cc(k);

[abserr,relerr] = get_err(dir_scattChebVec, dir.chebVecFromFourierVec(dir_scattFourierVec));
if verbose || abserr > 1e-9
    fprintf('Error in consistency between sols in Fourier and Cheb bases with dir BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dir_scattChebVec, dir_scattChebVec_axisymm);
if verbose || abserr > 1e-12
    fprintf('Error between dirBC and dirBC_axisymm sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

%% PML with dir BCs

pml         = pmlBC        (Nt,Nrs,Nr_out,pmlVs,coords,Rs,RPML,l,dl,d2l);
pml_axisymm = pmlBC_axisymm(Nt,Nrs,Nr_out,pmlVs,coords,Rs,RPML,l,dl,d2l);

%%
% In Fourier basis

pmlfull_scattFourierVec         = pml.solve_full_fc(k);
pmlfull_scattFourierVec_axisymm = pml_axisymm.solve_full_fc(k);

[abserr,relerr] = get_err(pmlfull_scattFourierVec, pmlfull_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm full sols in Fourier basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

pmlfull_scattChebVec         = pml.solve_full_cc(k);
pmlfull_scattChebVec_axisymm = pml_axisymm.solve_full_cc(k);

[abserr,relerr] = get_err(pmlfull_scattChebVec, pmlfull_scattChebVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm full sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

%% PML BC from Schur complement

%%
% In Fourier basis

pml_scattFourierVec         = pml.solve_fc(k);
pml_scattFourierVec_axisymm = pml_axisymm.solve_fc(k);
pml_scattValuesVec = pml.solve(k);

[abserr,relerr] = get_err(pml_scattFourierVec, pml_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm sols in Fourier basis: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = get_err(pml_scattValuesVec, true_scattValuesVec);
if verbose || abserr > 1e-7
    fprintf('Error between pmlBC and true sols in values basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

pml_scattChebVec         = pml.solve_cc(k);
pml_scattChebVec_axisymm = pml_axisymm.solve_cc(k);

[abserr,relerr] = get_err(pml_scattChebVec, dir.chebVecFromFourierVec(pml_scattFourierVec));
if verbose || abserr > 1e-8
    fprintf('Error in consistency between sols in Fourier and Cheb bases with dir BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(pml_scattChebVec, pml_scattChebVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between pmlBC and pmlBC_axisymm full sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

%% DtN BCs

dtn         = DtNBC        (Nt,Nrs,Vs,coords,Rs);
dtn_axisymm = DtNBC_axisymm(Nt,Nrs,Vs,coords,Rs);

%%
% In Fourier basis

dtn_scattFourierVec         = dtn.solve_fc(k);
dtn_scattFourierVec_axisymm = dtn_axisymm.solve_fc(k);
dtn_scattValuesVec = dtn.solve(k);

[abserr,relerr] = get_err(dtn_scattFourierVec, dtn_scattFourierVec_axisymm);
if verbose || abserr > 1e-11
    fprintf('Error between DtNBC and DtNBC_axisymm sols: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = get_err(dtn_scattValuesVec, true_scattValuesVec);
if verbose || abserr > 1e-8
    fprintf('Error between DtNBC and true sols in values basis: %4.2e, %4.2e\n', abserr, relerr);
end

%%
% In Cheb basis

dtn_scattChebVec         = dtn.solve_cc(k);
dtn_scattChebVec_axisymm = dtn_axisymm.solve_cc(k);

[abserr,relerr] = get_err(dtn_scattChebVec, dtn.chebVecFromFourierVec(dtn_scattFourierVec));
if verbose || abserr > 1e-9
    fprintf('Error in consistency between sols in Fourier and Cheb bases with DtN BCs: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dtn_scattChebVec, dtn_scattChebVec_axisymm);
if verbose || abserr > 1e-12
    fprintf('Error between DtNBC and DtNBC_axisymm sols in Cheb basis: %4.2e, %4.2e\n', abserr, relerr);
end

[abserr,relerr] = get_err(dtn.C(k), pml.C(k));
if verbose || abserr > 1e-5
    fprintf('Error between dtn.C and pml.C at k = %f: %4.2e, %4.2e\n', k, abserr, relerr);
end

warning on
disp('Warnings are back now');


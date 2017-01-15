function test_chebcoeff_picture
% Solving for scattered wave (as stack of fourier coefficients evaluated on
% mesh in r direction). The potential is piecewise constant and axisymmetric.

addpath(fileparts(pwd));
close all

%% Setup
%%
% choose piecewise constant, axisymmetric 2D potential and stipulate
% coordinates in which it is given
Rin = 0.88; Rout = 1.67;
coords = 'polar';
V0s = [0,50];
Vs = {@(r,t) 0*r + V0s(1); @(r,t) 0*r + V0s(2)};

%%
% pick frequency of incident wave
k = 2;
incfun = @(r,t) exp(1i*k*r.*cos(t));

%%
% choose discretization parameters
Ns = 20; 
Nrs = [25,25];
Rs = [Rin,Rout];

%% Discretize and solve

dtn = DtNBC(Ns,Nrs,Vs,coords,Rs);
RHS = dtn.RHSfromFun(Vs,incfun,coords);
scatt = dtn.T(k)\RHS;

%% Compare to analytically computed solution

truescatt = compute_truescatt(k,V0s,dtn);
%%
% compare
abserr = norm( scatt - truescatt ); relerr = abserr/norm( truescatt );
fprintf('Error for stack of fourier coeffs of scatt: %4.2e, %4.2e\n', abserr,relerr);

%% Cheb coeffs picture

RHS_cc = dtn.RHSfromFun_cc(Vs,incfun,coords);
scatt_cc = dtn.T_cc(k)\RHS_cc;

abserr = norm( dtn.fourierVecFromChebVec(scatt_cc) - scatt ); relerr = abserr/norm( scatt );
fprintf('Error in scatt in fourier basis: %4.2e, %4.2e\n', abserr,relerr);

abserr = norm( dtn.T(k)*dtn.fourierVecFromChebVec(scatt_cc) - RHS ); relerr = abserr/norm( RHS );
fprintf('Residual for Cheb basis solution w.r.t. Fourier basis equation: %4.2e, %4.2e\n', abserr, relerr);

abserr = norm( dtn.T_cc(k)*dtn.chebVecFromFourierVec(scatt) - RHS_cc ); relerr = abserr/norm( RHS_cc );
fprintf('Residual for Fourier basis solution w.r.t. Cheb basis equation: %4.2e, %4.2e\n', abserr, relerr);

%% Schur complement of cc -> cc operator
%%
% first make two problems from finer discretizations (not necessarily with
% same discretization points)
dNrs = 10;
dtn1 = DtNBC(dtn.Ns,dtn.Nrs+dNrs,dtn.Vs,coords,Rs);
eqtn1_cc = make_eqtn_cc(dtn1,dtn,k,Vs,incfun,coords);

dNrs = 15;
dtn2 = DtNBC(dtn.Ns,dtn.Nrs+dNrs,dtn.Vs,coords,Rs);
eqtn2_cc = make_eqtn_cc(dtn2,dtn,k,Vs,incfun,coords);

%%
% compare their Schur complements

abserr = norm( eqtn1_cc.Tk_comp - eqtn2_cc.Tk_comp );
relerr = abserr/norm( eqtn1_cc.Tk_comp );
fprintf('Error between schur comps of refined probs: %4.2e, %4.2e\n', abserr,relerr);

abserr = norm( scatt_cc - eqtn1_cc.scatt_comp ); relerr = abserr/norm( scatt_cc );
fprintf('Error between scattering solution and one obtained by solving\n');
fprintf('  Schur complemented version of first  refined system: %4.2e, %4.2e\n', abserr,relerr);

abserr = norm( scatt_cc - eqtn2_cc.scatt_comp ); relerr = abserr/norm( scatt_cc );
fprintf('Error between scattering solution and one obtained by solving\n');
fprintf('  Schur complemented version of second refined system: %4.2e, %4.2e\n', abserr,relerr);

end

function eqtn_cc = make_eqtn_cc(dtn,dtn_coarse,k,Vs,incfun,coords)
% set up a struct to keep track of different parts of equation and its
% relationship to coarse problem
    % set up the solve for the refined problem
    T_cc = dtn.T_cc(k);
    RHS_cc = dtn.RHSfromFun_cc(Vs,incfun,coords);
    scatt_cc = T_cc\RHS_cc;

    % subsetting 
    Isub = [1:dtn_coarse.Nrs(1), dtn.Nrs(1) + (1:dtn_coarse.Nrs(2))];
    I1 = repmat(              Isub',           1,dtn.Nt) + ...
         repmat((0:dtn.Nt-1)*dtn.Nr,length(Isub),     1);
    I1 = I1(:); I2 = 1:length(dtn.A); I2(I1) = [];
    sc = schurComp(I1,I2);

    % set up solve for Schur complemented problem
    T_cc_comp = sc.S(T_cc);
    RHS_cc_comp = sc.SRHS(T_cc,RHS_cc); 
    scatt_cc_comp = T_cc_comp\RHS_cc_comp;
    
    % make structure
    eqtn_cc = struct('Tk',T_cc,'RHS',RHS_cc,'scatt',scatt_cc,...
                     'Tk_comp',T_cc_comp,'RHS_comp',RHS_cc_comp,...
                     'scatt_comp',scatt_cc_comp);
end
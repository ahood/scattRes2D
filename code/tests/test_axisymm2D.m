function test_axisymm2D(verbose)

close all
addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end

disp('Temporarily turning warnings off');
warning off

do_plots = verbose; % whether we want to see plots or not
show_errors = verbose; % whether we want to see numerical values of errors
                       % regardless of whether they're small

% choose piecewise constant, axisymmetric 2D potential
Rin    = 1.5; 
Rout   = 2.5; 
V0_in  = 0; 
V0_out = 38;
coords = 'polar';
Vs = {@(r,t) 0*r + V0_in; @(r,t) 0*r + V0_out};
pot_params = struct('Vs',{Vs},'Rs',[Rin,Rout],'Vvals',[V0_in,V0_out]);

% choose discretization parameters and PML parameters
RPML = 8; 
Nr_out = 40;
a = 0.4;
l = @(t) a*(t-Rout).^3;
dl = @(t) a*3*(t-Rout).^2;
d2l = @(t) a*3*2*(t-Rout);
pmlVs = Vs; pmlVs{end+1} = @(r,t) 0*r; % isn't this wrong?
pml_params = struct('pmlVs',{pmlVs},'l',l,'dl',dl,'d2l',d2l,...
    'Nr_out',Nr_out,'RPML',RPML);

% Choosing discretization of potential support
Nt = 30;
Nrs = [20,20]; 

% Choosing energy of incident wave and number of fourier coefficients
z = 1.5;
maxn = 15;

% Scattering computation
do_dtn = true; do_pml = true; do_dir = true;
axisymm2D_scattering(pot_params,Nrs,Nt,pml_params,z,maxn,do_dtn,do_pml,do_dir,do_plots,show_errors);

warning on
disp('Warnings are back now');

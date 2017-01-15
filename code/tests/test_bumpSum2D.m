function test_bumpSum2D(verbose)

close all
addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end

% two bumps
d = 2; % center-center distance
r1 = d/2; r2 = d/4; % respective radii
c1 = (r1-r2)/2 - d/2;
c2 = (r1-r2)/2 + d/2; % respective centers (chosen so the smallest possible origin-centered disc encloses both)
V01 = 50; V02 = 30; % respective heights
V = @(x,y) bump( sqrt((x-real(c1)).^2 + (y-imag(c1)).^2)/r1 )*V01*exp(1) + ...
           bump( sqrt((x-real(c2)).^2 + (y-imag(c2)).^2)/r2 )*V02*exp(1);
Rs = c2+r2;       
Vs = {V}; coords = 'rect';

% choose discretization parameters and PML parameters
Nt = 16; Nrs = 40;
RPML = 6; Nr_out = 40;
a = 0.4;
l = @(t) a*(t-Rs).^3;
dl = @(t) a*3*(t-Rs).^2;
d2l = @(t) a*3*2*(t-Rs);
pmlVs = Vs; pmlVs{end+1} = @(x,y) 0*x;

% set up DtN, Dir and PML problems for comparison to one another
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dir = dirBC(Nt,Nrs,Vs,coords,Rs);
pml = pmlBC(Nt,Nrs,Nr_out,pmlVs,coords,Rs,RPML,l,dl,d2l);

% pick frequency of incident wave
z = 1.5; k = sqrt(z);

% compute scattering solutions for each type of BC
dtnscatt = dtn.solve(k);
pmlscatt = pml.solve(k);
dirscatt = dir.solve(k);

%% Checking it

% show potential as sanity check
if verbose
    figure, dtn.plotFun(dtn.Vs,'bumpSum2D potential',@real,coords);
end

[abserr,relerr] = L2err(dtn,dtnscatt,dir,dirscatt,[dtn.theta, 2*pi]);
if verbose || relerr > 1
    fprintf('L2 error between dtn and dir sols: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = L2err(dtn,dtnscatt,dtn,pmlscatt,[dtn.theta, 2*pi]);
if verbose || abserr > 1e-8
    fprintf('L2 error between dtn and pml sols: %4.2e, %4.2e\n', abserr, relerr);
end

end

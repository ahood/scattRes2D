function test_singleBump2D(verbose)

close all

addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end

% single bump function
Nt = 20; Nrs = 30;
Rs = 1; % radius of bump function
V0 = 50; V = @(x,y) bump(sqrt(x.^2+y.^2)/Rs)*V0*exp(1);
Vs = {V};
coords = 'rect';

% set up DtN, Dir and PML problems for comparison to one another
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dir = dirBC(Nt,Nrs,Vs,coords,Rs);

% use PML outside of Rout
RPML = 3*Rs;
a = 1;
l = @(t) a*(t - Rs).^3;
dl = @(t) a*3*(t - Rs).^2;
d2l = @(t) a*3*2*(t - Rs);
pmlVs = Vs; pmlVs{end+1} = @(x,y) 0*x;
pml = pmlBC(Nt,Nrs,Nrs,pmlVs,coords,Rs,RPML,l,dl,d2l);

% solve scattering problem in each case
z = 1.5; k = sqrt(z);

dtnscatt = dtn.solve(k);
pmlscatt = pml.solve(k);
dirscatt = dir.solve(k);

[abserr,relerr] = L2err(dtn,dtnscatt,dir,dirscatt,[dtn.theta,2*pi]);
if verbose || abserr > 1
    fprintf('L2 error between dtn and dir sols: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = L2err(dtn,dtnscatt,dtn,pmlscatt,[dtn.theta,2*pi]);
if verbose || abserr > 1e-8
    fprintf('L2 error between dtn and pml sols: %4.2e, %4.2e\n', abserr, relerr);
end

end

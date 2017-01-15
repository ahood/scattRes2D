function test_sqrBump2D(verbose)

close all

addpath(fileparts(pwd));
if nargin == 0, verbose = 0; end

Nt = 20; Nrs = 30; 
w = 1;     % width of bump function
p = 10;    % put centers of bump functions on lp ball
b = 1.1*w; %   of radius b (leaves some space in the middle)
Rs = sqrt(2)*b + w + 0.22; % radius at which DtN BC will be applied
V0 = 50; V = @(x,y) V0*squareBump(x,y,p,w,b)*exp(1); % so it's height V0
Vs = {V};
coords = 'rect';

% set up DtN, Dir and PML problems for comparison to one another
dtn = DtNBC(Nt,Nrs,Vs,coords,Rs);
dir = dirBC(Nt,Nrs,Vs,coords,Rs);

% use PML outside of Rout
RPML = 3*Rs;
a = 0.2;
l = @(t) a*(t - Rs).^3;
dl = @(t) a*3*(t - Rs).^2;
d2l = @(t) a*3*2*(t - Rs);
pmlVs = Vs; pmlVs{end+1} = @(x,y) 0*x;
pml = pmlBC(Nt,Nrs,40,pmlVs,coords,Rs,RPML,l,dl,d2l);

% solve scattering problem in each case
z = 1.5; k = sqrt(z); 

dtnscatt = dtn.solve(k);
pmlscatt = pml.solve(k);
dirscatt = dir.solve(k);

[abserr,relerr] = L2err(dtn,dtnscatt,dir,dirscatt,[dtn.theta,2*pi]);
if verbose || relerr > 1
    fprintf('L2 error between dtn and dir sols: %4.2e, %4.2e\n', abserr, relerr);
end
[abserr,relerr] = L2err(dtn,dtnscatt,dtn,pmlscatt,[dtn.theta,2*pi]);
if verbose || abserr > 1e-9
    fprintf('L2 error between dtn and pml sols: %4.2e, %4.2e\n', abserr, relerr);
end

end

function Vxy = squareBump(x,y,p,w,b)
    t = atan2(y,x); r = sqrt(x.^2 + y.^2);
    c = b./(cos(t).^p + sin(t).^p).^(1/p);
    Vxy = bump((r-c)/w);
    Vxy( isnan(Vxy) ) = 0;
end

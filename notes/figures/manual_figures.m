close all
addpath('../../code');
ts = linspace(0,2*pi,100);

% simple mesh
Nt = 10; Nr = 5; R = 1;
s = scattResComp2d(Nt,Nr,R);

figure, hold on
plot(s.xs,s.ys,'*');
plot(exp(1i*ts),'g');
axis([-1.5 1.5 -1.5 1.5]);
axis square
box on
print -depsc simplemesh.eps

% complicated mesh
Nt = 16; Nrs = [3,4,5]; Rs = [0.5,1,2];
s = scattResComp2d(Nt,Nrs,Rs);

figure, hold on
plot(s.xs,s.ys,'*');
plot(0.5*exp(1i*ts),'g');
plot(    exp(1i*ts),'g');
plot(  2*exp(1i*ts),'g');
axis([-2.5 2.5 -2.5 2.5]);
axis square
box on
print -depsc complicatedmesh.eps

% axisymmetric, piecewise constant ring
Nt = 30; Nrs = [10,10]; Rs = [1,2];
s = scattResComp2d(Nt,Nrs,Rs);

V = @(r,t) (1 < r) & (r < 2);

figure, s.plotFun(V,[],[],'polar');
print -depsc ringpotential.eps

figure, plot(s.xs,s.ys,'*');
axis([-2.5 2.5 -2.5 2.5]);
axis square
box on
print -depsc ringpotentialmesh.eps

% three bumps
Nt = 50; Nrs = 30; Rs = 3;
s = scattResComp2d(Nt,Nrs,Rs);

Rcenters = 1.4;
sigma = 1/3;
c = Rcenters*exp(2i*pi*(1:3)/3);
G = @(x,y,cj) exp(-abs(x+1i*y-cj).^2/2/sigma^2);
V = @(x,y) (G(x,y,c(1)) + G(x,y,c(2)) + G(x,y,c(3)));

figure, s.plotFun(V,[],[],'rect');
print -depsc threebumppotential.eps

figure, plot(s.xs,s.ys,'*');
axis([-3.5 3.5 -3.5 3.5]);
axis square
box on
print -depsc threebumppotentialmesh.eps

% PML
sigma = @(t) 0.5*(t-Rs(end)).*(t >= Rs(end));
dsigma = @(t) 0.5.*(t >= Rs(end));

l = @(t) 0.25*(t-Rs(end)).^2 .*(t >= Rs(end));
dl = @(t) sigma(t);
d2l = @(t) dsigma(t);
r2s = @(r) r + 1i*l(r);

pml = pmlBC(Nt,Nrs,7,{V,@(x,y) 0*x},'rect',Rs,5,l,dl,d2l);

figure, hold all
plot(pml.r, 0*pml.r, '-*', 'displayname', 'r');
plot(pml.s, '-*', 'displayname', 's(r)');
xlabel('real')
ylabel('imag')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1:3) pos(4)*3/4]);
legend show
legend('location','best')
box on

print -depsc pml.eps


close all
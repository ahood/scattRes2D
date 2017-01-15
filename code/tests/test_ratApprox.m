function test_ratApprox(verbose)
% experimenting with getting rational approximation to function with poles

close all
addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end

%% Example 1: analytic scalar-valued function

% function
n = 5;
f = @(z) besselh(n,1,sqrt(z)*5);
ell = ellipse(10,pi/4,5,3,[],[]);
r = rect(5,15,-5,5,100,100);
t = linspace(0,1,100);
ratf = ratApprox(ell,f,t);
ratf.show_error(r);
title('test\_ratApprox: Error in rat approx to besselh(5,1,sqrt(z)*5)');

%% Example 2: function with a few poles

% function
p1 = 10; p2 = 11; p3 = 10+1i;
f = @(z) sin(z)./(z-p1)./(z-p2)./(z-p3);

poles = [p1, p2, p3];
rs = 0*(1:3) + .05;
c = circ(10,2,[],[]);
r = rect(8,12,-2,2,100,100);
t = linspace(0,1,200);

ratf = ratApproxWithPoles(c,f,poles,rs,t);
ratf.show_error(r);
title('test\_ratApprox: Error in rat approx to sin(z)/(z-10)/(z-11)/(z-(10+1i))');

%% Example 3: one component of DtN map

% define DtN map function
n = 8;
Rad = 1.75;
f = @(z) DtNBC.DtNcoeffs(n,sqrt(z),Rad);
r = rect(1,50,-50,0,100,100);
r.log10contour(@(z) abs(f(z))); 
title('test\_ratApprox: Log10 magnitude of DtNcoeffs(8,sqrt(z),1.75)');

g  = @(z) scattResComp2d.H(n,sqrt(z),Rad);
dg = @(z) scattResComp2d.dH(n,sqrt(z),Rad)*Rad/2/z;
z0 = 8.86-10.65i;
for jj = 1:5
    z0 = z0 - g(z0)/dg(z0);
end

poles = [z0];
rs = 0*poles + 0.05;
c = circ(15-10i,10,[],[]);
r = rect(5,25,-20,0,100,100);
t = linspace(0,1,1000);

ratf = ratApproxWithPoles(c,f,poles,rs,t);
ratf.show_error(r);
hold on, plot(ratf.poles,'r*');
title('test\_ratApprox: Error in rat approx to DtNcoeffs(8,sqrt(z),1.75) (poles shown)');

%% Example 4: Many poles over large region

r = linspace(1,10,10);
t = linspace(0,2*pi,10); t = t(1:end-1);
[rr,tt] = meshgrid(r,t);
poles = rr.*cos(tt) + 1i*rr.*sin(tt); poles = poles(:).';

function fevalz = feval(z)
    fevalz = length(poles)*sin(z);
    for jj = 1:length(poles)
        fevalz = fevalz./(z-poles(jj));
    end
end

rs = 0*poles + .05;
c = circ(0,11,[],[]);
rec = rect(-12,12,-12,12,100,100);
t = linspace(0,1,200);

ratf = ratApproxWithPoles(c,@feval,poles,rs,t);
ratf.show_error(rec);
hold on, plot(ratf.poles,'r*');
title('test\_ratApprox: Error in rat approx to function with lots of poles (poles shown)');

end

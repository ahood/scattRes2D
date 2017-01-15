function test_L2norm_stuff

close all
addpath(fileparts(pwd));

%% test that a function on two different meshes has zero L2 norm
f = @(r,theta) exp(1i*r.*cos(theta));

% radius of disc
R = 3.9;

% first mesh and values
Nt1 = 20; %240;
Nr1 = 950; %0;
s1 = scattResComp2d_sparse(Nt1,Nr1,R);
f_values1 = f(s1.rs, s1.ts);

% second mesh and values
Nt2 = 10;
Nr2 = 245;
s2 = scattResComp2d_sparse(Nt2,Nr2,R);
f_values2 = f(s2.rs, s2.ts);

% theta mesh for comparison
other_theta = linspace(0,2*pi,11);

% L2 distance between the two function representations
[abserr,relerr,~] = L2err(s1,f_values1,s2,f_values2,other_theta);
fprintf('Error between same function on different meshes: %4.2e, %4.2e\n', abserr, relerr);

%% test L2 norm of f(r,theta) = r*sin(theta). 
R = 5;
theta = linspace(0,2*pi,21); % odd because mesh without 0 is even size
f_ontheta = cell(1,length(theta));
for ii = 1:length(theta)
    f_ontheta{ii} = chebfun(@(r) r*sin(theta(ii)), [-R,R]);
end
output = L2normOnDisc(f_ontheta,theta);
expected = sqrt(pi)*R^2/2;
abserr = abs(output-expected);
relerr = abserr/abs(expected);
fprintf('Error in L2-norm of rsin(theta) on B(0,R): %4.2e, %4.2e\n', abserr, relerr);

%% test L2 norm of f(r,theta) = r*exp(1i*theta)
R = 3;
theta = linspace(0,2*pi,15); % odd because mesh without 0 is even size
f_ontheta = cell(1,length(theta));
for ii = 1:length(theta)
    f_ontheta{ii} = chebfun(@(r) r*exp(1i*theta(ii)), [-R,R]);
end
output = L2normOnDisc(f_ontheta,theta);
expected = sqrt(pi/2)*R^2;
abserr = abs(output-expected);
relerr = abserr/abs(expected);
fprintf('Error in L2-norm of rexp(i theta) on B(0,R): %4.2e, %4.2e\n', abserr, relerr);

%% test L2 norm of f(r,theta) = r*exp(1i*theta)
R = 3;
theta = linspace(0,2*pi,3); % odd because mesh without 0 is even size
f_ontheta = cell(1,length(theta));
for ii = 1:length(theta)
    f_ontheta{ii} = chebfun(@(r) exp(1i*r*cos(theta(ii))), [-R,R]);
end
output = L2normOnDisc(f_ontheta,theta);
expected = sqrt(pi)*R;
abserr = abs(output-expected);
relerr = abserr/abs(expected);
fprintf('Error in L2-norm of exp(i x) on B(0,R): %4.2e, %4.2e\n', abserr, relerr);

%% test mapping function of theta onto finer theta mesh
g = @(t) 3*cos(t).^2;
theta = linspace(0,2*pi,6); theta = theta(1:end-1);
g_values = g(theta);

Nt = length(theta);
c = fft(g_values');

other_theta = linspace(0,2*pi,21)'; % don't need other_theta(1) = 0 or anything
g_values_other_theta = mapToOtherThetaMesh(g_values,1,other_theta);
abserr = norm(g_values_other_theta - g(other_theta));
relerr = abserr/norm(g(other_theta));
fprintf('Error in function of theta vals at finer mesh: %4.2e, %4.2e\n', abserr, relerr);

%% test mapping function of r and theta onto finer theta mesh
Nt = 40; Nrs = [20]; Rs = [3]; 
s = scattResComp2d_sparse(Nt,Nrs,Rs);
Nt_comp = 46;
s_comp = scattResComp2d_sparse(Nt_comp,Nrs,Rs);
other_theta = s_comp.theta;

f = @(r,t) exp(1i*r.*cos(t));
f_values = f(s.rs,s.ts);
f_values_other_theta = f(s_comp.rs,s_comp.ts);
abserr = norm(f_values_other_theta - mapToOtherThetaMesh(f_values,s.Nr,other_theta));
relerr = abserr/norm(f_values_other_theta);
fprintf('Error in function of r and theta vals at finer mesh: %4.2e, %4.2e\n', abserr, relerr);

% figure, s_comp.plotValuesVec(mapToOtherThetaMesh(f_values,s.Nr,other_theta),'mapped',@real);
% figure, s_comp.plotValuesVec(f_values_other_theta,'true',@real);

%% test L2 distance between two functions on different theta meshes
f = @(r,t) 3*r.^2.*cos(t);
sf = scattResComp2d_sparse(30,[20],[3]);
f_values = f(sf.rs,sf.ts);

g = @(r,t) r.*sin(t) + f(r,t);
sg = scattResComp2d_sparse(36,[20],[3]);
g_values = g(sg.rs,sg.ts);

other_theta = linspace(0,2*pi,61);
f_values_other_theta = mapToOtherThetaMesh(f_values,sf.Nr,other_theta);
g_values_other_theta = mapToOtherThetaMesh(g_values,sf.Nr,other_theta);

h_values_other_theta = f_values_other_theta - g_values_other_theta;
h_rect = reshape(h_values_other_theta,sf.Nr,[]);
h_radial_chebfuns = cell(1,length(other_theta));
for ii = 1:length(other_theta)
    h_radial_chebfuns{ii} = chebfun([flipud(h_rect(:,ii)); h_rect(:,ii)], [-sf.Rs(end), sf.Rs(end)]);
end
output = L2normOnDisc(h_radial_chebfuns, other_theta);
expected = sqrt(pi)*sf.Rs(end)^2/2;
abserr = abs(output-expected);
relerr = abserr/abs(expected);
fprintf('Error for L2-norm difference of two functions on different theta meshes of B(0,R): %4.2e, %4.2e\n', abserr, relerr);
[abserr2,relerr2] = L2err(sf,f_values,sg,g_values,other_theta);
fprintf('Error in L2err comp: %4.2e\n', abs(abserr2 - expected));

%% test L2 distance between two functions on different theta and r meshes
f = @(r,t) 3*r.^2.*cos(t);
sf = scattResComp2d_sparse(30,[20],[3]);
f_values = f(sf.rs,sf.ts);

g = @(r,t) r.*sin(t) + f(r,t);
sg = scattResComp2d_sparse(36,[25],[3]);
g_values = g(sg.rs,sg.ts);

other_theta = linspace(0,2*pi,61);
f_values_other_theta = mapToOtherThetaMesh(f_values,sf.Nr,other_theta);
g_values_other_theta = mapToOtherThetaMesh(g_values,sg.Nr,other_theta);

f_rect = reshape(f_values_other_theta,sf.Nr,[]);
g_rect = reshape(g_values_other_theta,sg.Nr,[]);

h_radial_chebfuns = cell(1,length(other_theta));
for ii = 1:length(other_theta)
    f_chebfun = chebfun([flipud(f_rect(:,ii)); f_rect(:,ii)], [-sf.Rs(end), sf.Rs(end)]);
    g_chebfun = chebfun([flipud(g_rect(:,ii)); g_rect(:,ii)], [-sf.Rs(end), sf.Rs(end)]);
    h_radial_chebfuns{ii} = f_chebfun - g_chebfun;
end
output = L2normOnDisc(h_radial_chebfuns, other_theta);
expected = sqrt(pi)*sf.Rs(end)^2/2;
abserr = abs(output-expected);
relerr = abserr/abs(expected);
fprintf('Error for L2-norm difference of two functions on different meshes of B(0,R): %4.2e, %4.2e\n', abserr, relerr);
[abserr2,relerr2] = L2err(sf,f_values,sg,g_values,other_theta);
fprintf('Error in L2err comp: %4.2e\n', abs(abserr2 - expected));

% finally, very important, check that when we represent the same function
% on two different meshes that the L2 error between the two representations
% is zero
f = @(r,t) 3*r.^2.*cos(t);

% reference solution from rep_lin_paper_sparse
Nr = 250;
Nt = 240;
R = 3.9; 
sf_ref = scattResComp2d_sparse(Nt,Nr,R);
f_values_ref = f(sf_ref.rs, sf_ref.ts);

% mesh for comparison
other_theta = linspace(0,2*pi,Nt+1); 

% compare with representations on various coarser meshes
Nrs = 30:50:230;
Nts = 50:50:250; % one line for each Nt
open('L2_err_disc_dependence.fig');
% figure, hold all
% for Nt = Nts
%     fprintf('Working on Nt = %d\n', Nt);
%     L2errs = 0*Nrs;
%     for ii = 1:length(Nrs)
%         Nr = Nrs(ii);
%         fprintf('--Working on Nr = %d\n', Nr);
%         tic;
%         sf = scattResComp2d_sparse(Nt,Nr,R);
%         fprintf('  Took %.2f s to make scattResComp2d_sparse instance\n', toc);
%         f_values = f(sf.rs, sf.ts);
%         tic;
%         L2errs(ii) = L2err(sf_ref, f_values_ref, sf, f_values, other_theta);
%         fprintf('  Took %.2f s to compute L2 error with reference\n', toc);
%     end
%     plot(Nrs, log10(L2errs), 'displayname', ['Nt=' num2str(Nt)]);
% end
% xlabel(Nr);
% ylabel('log10 error');
% title(['Log 10 absolute error with same function evaluated on Nr=' num2str(Nr) ', Nt=' num2str(Nt)])

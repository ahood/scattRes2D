function laplacian_construction

close all 
addpath(fileparts(pwd));

%% Setup for half Chebyshev mesh in r direction

% choose size of meshes in r and theta directions
Nt = 4; % must be even
Nr = 3;

% theta mesh and second deriv matrix
dtheta = 2*pi/Nt;
theta = (dtheta:dtheta:2*pi)';
e = ones(Nt,1);
D2t = toeplitz([-pi^2/(3*dtheta^2)-1/6 ...
                 .5*(-1).^(2:Nt)./sin(dtheta*(1:Nt-1)/2).^2]);
It = eye(Nt);             
Ithalf = eye(Nt/2);

% r mesh and first and second deriv matrices
% [D2r,Dr,r] = radialCheb(Nr,0,1,0);
[D,x] = cheb(2*Nr-1); % vector length 2*Nr returned
r = x(Nr+1:end);
D2 = D*D;
% blocks we need
D3  = D2(Nr+1:end,1:Nr);
D3t = D2(Nr+1:end,Nr:-1:1);
D4  = D2(Nr+1:end,Nr+1:end);
E3  =  D(Nr+1:end,1:Nr);
E3t =  D(Nr+1:end,Nr:-1:1);
E4  =  D(Nr+1:end,Nr+1:end);

R = diag(1./r); % needed later
Ir = eye(Nr);

% discrete laplace operator pieces
J = [0*Ithalf, Ithalf; Ithalf, 0*Ithalf];
Dr  = kron(It,E4) + kron(J,E3t);
Drr = kron(It,D4) + kron(J,D3t);
Dtt = kron(D2t,Ir);
L = kron(It,D4 + R*E4) + kron(J,D3t + R*E3t) + kron(D2t,R*R);

% check that pieces add to L
abserr = norm( Drr + kron(It,R)*Dr + kron(It,R*R)*Dtt - L );
relerr = abserr/norm(L);
fprintf('Error between laplacian and sum of pieces: %4.2e, %4.2e\n', ...
    abserr, relerr);

% augmented meshes
r_aug = repmat(r,Nt,1);
theta_aug = ones(size(r))*theta.'; theta_aug = theta_aug(:);

%% Test for half Chebyshev mesh in r direction

% define a function and its laplacian, and compare to discretized
% approximation of laplacian

% u = x + y
% u    = @(r,t) r.*cos(t) + r.*sin(t);
% u_r  = @(r,t) cos(t) + sin(t);
% u_rr = @(r,t) 0*r;
% u_t  = @(r,t) r.*-sin(t) + r.*cos(t);
% u_tt = @(r,t) r.*-cos(t) + r.*-sin(t);

% u = x^2 + y^2
u    = @(r,t) (r.*cos(t)).^2 + (r.*sin(t)).^2;
u_r  = @(r,t) 2*(r.*cos(t)).*cos(t) + ...
              2*(r.*sin(t)).*sin(t);
u_rr = @(r,t) 2*cos(t).^2 + 2*sin(t).^2;
u_t  = @(r,t) 0*r; % u is a function of r
u_tt = @(r,t) 0*r;          

u_true = u(r_aug,theta_aug);     

% check u_r
abserr = norm( Dr*u_true - u_r(r_aug,theta_aug) );
relerr = abserr/norm(u_r(r_aug,theta_aug));
fprintf('Error in u_r action: %4.2e, %4.2e\n', abserr,relerr);

% check u_rr
abserr = norm( Drr*u_true - u_rr(r_aug,theta_aug) );
relerr = abserr/norm(u_rr(r_aug,theta_aug));
fprintf('Error in u_rr action: %4.2e, %4.2e\n', abserr,relerr);

% check u_tt
abserr = norm( Dtt*u_true - u_tt(r_aug,theta_aug) );
relerr = abserr/norm(u_tt(r_aug,theta_aug));
fprintf('Error in u_tt action: %4.2e, %4.2e\n', abserr, relerr);

% check laplacian
lap_true = u_rr(r_aug,theta_aug) + u_r(r_aug,theta_aug)./r_aug + ...
           u_tt(r_aug,theta_aug)./r_aug.^2;
abserr = norm(L*u_true - lap_true);
relerr = abserr/norm(lap_true);
fprintf('Error in laplacian: %4.2e, %4.2e\n', abserr, relerr);

%% Setup for piecewise meshes

% choose size of meshes in r and theta directions
Nt = 4; % must be even
Nrs = [3,5,7]; % three meshes
Rs = [0.75, 1.2, 2.4]; % on intervals with these boundaries

% theta mesh and second deriv matrix
dtheta = 2*pi/Nt;
theta = (dtheta:dtheta:2*pi)';
e = ones(Nt,1);
D2t = toeplitz([-pi^2/(3*dtheta^2)-1/6 ...
                 .5*(-1).^(2:Nt)./sin(dtheta*(1:Nt-1)/2).^2]);
It = eye(Nt);             
Ithalf = eye(Nt/2);

%%
% Half Chebyshev mesh on [0,Rs(1)] 

% r mesh and first and second deriv matrices
[D,x] = cheb(2*Nrs(1)-1,-Rs(1),Rs(1)); % vector length 2*Nr returned
r1 = x(Nrs(1)+1:end);
D2 = D*D;
% blocks we need
D3  = D2(Nrs(1)+1:end,1:Nrs(1));
D3t = D2(Nrs(1)+1:end,Nrs(1):-1:1);
D4  = D2(Nrs(1)+1:end,Nrs(1)+1:end);
E3  =  D(Nrs(1)+1:end,1:Nrs(1));
E3t =  D(Nrs(1)+1:end,Nrs(1):-1:1);
E4  =  D(Nrs(1)+1:end,Nrs(1)+1:end);
R1 = diag(1./r1);

%%
% Regular Chebyshev mesh on [Rs(1),Rs(2)]
[D,x] = cheb(Nrs(2)-1,Rs(1),Rs(2));
r2 = x;
Dr2 = D;
Drr2 = D^2;
R2 = diag(1./r2);

%%
% Regular Chebyshev mesh on [Rs(2),Rs(3)]
[D,x] = cheb(Nrs(3)-1,Rs(2),Rs(3));
r3 = x;
Dr3 = D;
Drr3 = D^2;
R3 = diag(1./r3);

%%
% Combining
r = [r1; r2; r3];
R = diag(1./r);

%%
% discrete laplace operator pieces
J = [0*Ithalf, Ithalf; Ithalf, 0*Ithalf];
Dr  = kron(It,blkdiag( E4,  Dr2,  Dr3)) ...
    + kron( J,blkdiag(E3t,0*Dr2,0*Dr3));
Drr = kron(It,blkdiag( D4,  Drr2,  Drr3)) ...
    + kron( J,blkdiag(D3t,0*Drr2,0*Drr3));
L = kron(It,blkdiag( D4 + R1* E4,Drr2 + R2*Dr2,Drr3 + R3*Dr3)) ...
  + kron( J,blkdiag(D3t + R1*E3t,        0*Dr2,        0*Dr3)) ...
  + kron(D2t,R*R);

% augmented meshes
r_aug = repmat(r,Nt,1);
theta_aug = ones(size(r))*theta.'; theta_aug = theta_aug(:);

%% Test for piecewise meshes in r direction

fprintf('Running piecewise mesh tests...\n');

% define a function and its laplacian, and compare to discretized
% approximation of laplacian

% u = x + y
% u    = @(r,t) r.*cos(t) + r.*sin(t);
% u_r  = @(r,t) cos(t) + sin(t);
% u_rr = @(r,t) 0*r;
% u_t  = @(r,t) r.*-sin(t) + r.*cos(t);
% u_tt = @(r,t) r.*-cos(t) + r.*-sin(t);

% u = x^2 + y^2
u    = @(r,t) (r.*cos(t)).^2 + (r.*sin(t)).^2;
u_r  = @(r,t) 2*(r.*cos(t)).*cos(t) + ...
              2*(r.*sin(t)).*sin(t);
u_rr = @(r,t) 2*cos(t).^2 + 2*sin(t).^2;
u_t  = @(r,t) 0*r; % u is a function of r
u_tt = @(r,t) 0*r;          

u_true = u(r_aug,theta_aug);     

% check u_r
abserr = norm( Dr*u_true - u_r(r_aug,theta_aug) );
relerr = abserr/norm(u_r(r_aug,theta_aug));
fprintf('Error in u_r action: %4.2e, %4.2e\n', abserr,relerr);

% check u_rr
abserr = norm( Drr*u_true - u_rr(r_aug,theta_aug) );
relerr = abserr/norm(u_rr(r_aug,theta_aug));
fprintf('Error in u_rr action: %4.2e, %4.2e\n', abserr,relerr);

% check laplacian
lap_true = u_rr(r_aug,theta_aug) + u_r(r_aug,theta_aug)./r_aug + ...
           u_tt(r_aug,theta_aug)./r_aug.^2;
abserr = norm(L*u_true - lap_true);
relerr = abserr/norm(lap_true);
fprintf('Error in laplacian: %4.2e, %4.2e\n', abserr, relerr);


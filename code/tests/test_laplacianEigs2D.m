function test_laplacianEigs2D(verbose)
% Eigs and modes of Laplacian checked against Spectral Methods in MATLAB

close all
addpath(fileparts(pwd));

if nargin == 0, verbose = 0; end

% Nt = 20; Nr = 20;
Nt = 20; Nrs = [10, 15, 10]; Rs = [0.25,0.8,1];
my = dirBC(Nt,Nrs,{@(x,y) 0*x, @(x,y) 0*x, @(x,y) 0*x},'rect',Rs);
my.eig_comp();

indx = [1 3 6 10];
theirs = [1; 1.5933405057; 2.2954172674; 2.9172954551];
l = my.evals; V = my.evecs;
l = l(indx); V = V(:,indx);
mine = sqrt(l/l(1));

abserr = norm(mine-theirs);
relerr = abserr/norm(theirs);

if verbose || abserr > 1e-9
    fprintf('Error in computed eigs: %4.2e, %4.2e\n', abserr, relerr);
    
    z = exp(1i*pi*(-100:100)/100);
    [ay,ax] = meshgrid([.58 .1],[.1 .5]); clf

    for ii = 1:length(indx)
        subplot('position',[ax(ii) ay(ii) .4 .4])
        plot(z), axis(1.05*[-1 1 -1 1 -1 1]), axis off, hold on
        u = my.plotValuesVec(V(:,ii),'',@real);
        view(0,20), colormap([0 0 0]), axis square
        contour3(my.xx,my.yy,u-1,[-1 -1])
        plot3(real(z),imag(z),-abs(z))
        text(-.8,4,['Mode ' int2str(indx(ii))],'fontsize',9)
        text(-.8,3.5, ['theirs = ' num2str(theirs(ii),'%16.10f')],'fontsize',9)
        text(-.8,3  , ['mine   = ' num2str(mine(ii),  '%16.10f')],'fontsize',9)
    end
    figure, imshow(imread('SMM_laplacian_modes.jpg'));
end
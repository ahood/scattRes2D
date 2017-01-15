function [abserrs,relerrs] = axisymm2D_scattering(pot_params,Nrs,Nt,pml_params,z,maxn,do_dtn,do_pml,do_dir,do_plots,show_errors)
% Do the scattering computation with the specified boundary conditions.

if nargin < 11, show_errors = true; end
if nargin < 10, do_plots = false; end

% close all
addpath(fileparts(pwd));

k = sqrt(z);
coords = 'polar';
s = scattResComp2d(Nt,Nrs,pot_params.Rs); % need for pml soln plot

% analytical solution (approximate, since comes from truncating Fourier
truescatt = compute_truescatt(k,pot_params.Vvals,s,maxn);

% compute scattering solutions for each type of BC
abserrs = [];
relerrs = [];
other_theta = linspace(0,2*pi,s.Nt+1);
if do_dtn
    dtn = DtNBC(Nt,Nrs,pot_params.Vs,coords,pot_params.Rs);
    dtnscatt = dtn.solve(k);     dtnrect = reshape(dtnscatt,s.Nr,[]);
    [abserr,relerr] = L2err(s, truescatt, dtn, dtnscatt, other_theta);
    abserrs = [abserrs abserr]; relerrs = [relerrs relerr];
    if show_errors || abserr > 1e-9
        fprintf('For Nt = %d, Nr = %d...\n', s.Nt, s.Nr);
        fprintf('L2 error between true and dtn valuesVec: %4.2e, %4.2e\n', abserr,relerr);
    end
end
if do_pml
    pml = pmlBC(Nt,Nrs,pml_params.Nr_out,pml_params.pmlVs,...
        coords,pot_params.Rs,pml_params.RPML,pml_params.l,...
        pml_params.dl,pml_params.d2l); 
    pmlscatt = pml.solve(k);    pmlrect = reshape(pmlscatt,s.Nr,[]);
    [abserr,relerr] = L2err(s, truescatt, s, pmlscatt, other_theta);
    abserrs = [abserrs abserr]; relerrs = [relerrs relerr];
    if show_errors || abserr > 1e-7
        fprintf('L2 error between true and pml valuesVec: %4.2e, %4.2e\n', abserr,relerr);
    end
end
if do_dir
    dir = dirBC(Nt,Nrs,pot_params.Vs,coords,pot_params.Rs);
    dirscatt = dir.solve(k);    dirrect = reshape(dirscatt,s.Nr,[]);
    [abserr,relerr] = L2err(s, truescatt, dir, dirscatt, other_theta);
    abserrs = [abserrs abserr]; relerrs = [relerrs relerr];
    if show_errors || abserr > 1
        fprintf('L2 error between true and dir valuesVec: %4.2e, %4.2e\n', abserr,relerr);
    end
end

if do_dtn && do_pml
    [abserr,relerr] = L2err(dtn, dtnscatt, dtn, pmlscatt, other_theta);
    if show_errors || abserr > 1e-8
        fprintf('L2 error between dtn and pml valuesVec: %4.2e, %4.2e\n', abserr, relerr);
    end
end

%%
% Visualize absolute errors

if do_plots 
    if do_dir
        figure, s.plotValuesVec(log10(abs(dirscatt-truescatt)),'dir abserr');
    end
    if do_pml
        figure, s.plotValuesVec(log10(abs(pmlscatt-truescatt)),'pml abserr');
    end
    if do_dtn
        figure, s.plotValuesVec(log10(abs(dtnscatt-truescatt)),'dtn abserr');
    end
end

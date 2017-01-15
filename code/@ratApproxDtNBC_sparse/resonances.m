function [resvecs,resvals] = resonances(obj,n,E0,verbose)
    if nargin < 4, verbose = 0; end
    if nargin < 3 || isempty(E0), E0 = obj.ell.c; end
    if nargin < 2 || isempty(n), n = 40; end
    matsize = obj.dtns.Nr*obj.dtns.Nt + length(obj.A22diag);
    itmax = matsize;
    k0 = obj.mysqrt(E0); % pick the right branch
    T_k0 = @(x) obj.apply_Tfull(x,k0);
    M_k0 = @(x) obj.apply_pc_full(x,k0); % default preconditioner
    function y = tmpFun(x)
        [y,flag] = gmres(T_k0,x,[],1e-6,itmax,M_k0);
    end
    if verbose
        T_k0_inv = @(x) gmres(T_k0, x,[],1e-6,itmax,M_k0);
    else
        T_k0_inv = @tmpFun;
    end
    eigs_fun = @(x) T_k0_inv(obj.apply_Bfull(x));

    % compute 
    opts.isreal = false;
    if verbose
        fprintf('--Computing %d eigs closest to %s\n', n, num2str(k0));
        tic;
        [resvecs,resvals] = eigs(eigs_fun, matsize, n, E0,opts);
        fprintf('  --Took %.1f s\n', toc);
    else
        [resvecs,resvals] = eigs(eigs_fun, matsize, n, E0,opts);
    end
    resvals = diag(resvals);
end

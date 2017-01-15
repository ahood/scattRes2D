classdef ratApproxDtNBC_sparse < handle
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with approximate DtN boundary conditions,
    % where approximations are meant to be good in region
    properties
        br = 1; % branch of sqrt (1 or 2)
        dtns % DtNBC_sparse problem (I think this will just be a reference since subclass of handle, so no memory issue?)
        ell % ellipse in which the rational approximation is supposed to be pretty good
        mysqrt % either +sqrt or -sqrt
        ratf = []; % list of rational approx objects for individual components of DtN map
        ratf_lens = []; % size of each block
        ratf_ends % last row of each block
        poles = {};
        sc % schur complement object for map from scattered wave to -V*incident
        A22diag

        blk_precond_k % the value of k last used for preconditioner
        M11_blks % a cell array where j-th element is struct with L,U decomp of j-th block of block diag preconditioner M  
        blk_fourier_precond_k
        fourierS_blks % a cell array j-th elt is struct with L,U decomp of j-th block
        sc_blk_fourier_precond_k
        sc_fourierM_blks
    end
    methods
        function obj = ratApproxDtNBC_sparse(dtns,ell,br,N)
        % V is a potential function handle. V must be smooth on the R
        % disk and V = V(x,y), not V = V(r,theta).
            if nargin < 4, N = 500; end   
            obj.br = br; % 1 or 2
            obj.dtns = dtns;
            obj.ell = ell;
            
            % use given branch of square root function
            if br == 1, mysqrt = @(z) sqrt(z);
            elseif br == 2, mysqrt = @(z) -sqrt(z);
            end
            obj.mysqrt = mysqrt;
            
            % do first pass to get discretizations and sizes
            for jj = 1:dtns.Nt
                n = dtns.Ns(jj);
                fn = @(z) DtNBC.DtNcoeffs(n,mysqrt(z),dtns.r(end));
                poles{jj} = ratApproxDtNBC.get_poles_in_range(n,ell,dtns.r(end),br);
                t = linspace(0,1,N); % wow, need this to be better
                rs = 0*poles{jj} + 0.05; % wow, also so dumb
                ratfn = ratApproxWithPoles(ell,fn,poles{jj},rs,t);
                ratf(jj) = ratfn;
                ratf_lens(jj) = length(ratfn.w);
            end
            
            l = 0;
            for jj = 1:dtns.Nt
                cols = l + (1:ratf_lens(jj));
                A22diag( cols) = ratfn.z;
                l = l + ratf_lens(jj);
            end
            
            obj.ratf = ratf;
            obj.poles = poles;
            obj.ratf_lens = ratf_lens;
            obj.ratf_ends = cumsum(ratf_lens);

            obj.A22diag = A22diag(:);

            obj.sc = schurComp(1:(dtns.Nt*dtns.Nr), (dtns.Nt*dtns.Nr) + (1:length(A22diag)) );
            
            obj.blk_precond_k = Inf;
            obj.M11_blks = cell(1,obj.dtns.Nt);
            
            obj.blk_fourier_precond_k = Inf;
            obj.fourierS_blks = cell(1,obj.dtns.Nt);
            
            obj.sc_blk_fourier_precond_k = Inf;
            obj.sc_fourierM_blks = cell(1,obj.dtns.Nt);
        end
        
        % compute n resonances closest to E0
        [resvecs,resvals] = resonances(obj,n,E0,verbose)

        % solve the scattering problem on disc
        scattValuesVec = solve(obj,k,incfun,precond,verbose)
                
        % look at poles
        show_poles(obj,r)
        
        % don't form, apply 
        x = apply_T(obj,x0,k)
        x = apply_Tfull(obj,x0,k)
        x = apply_Bfull(obj,x0)
        x = apply_pc_blkdiag(obj,x0,k) % for Tfull
        x = apply_pc_fourier_Vav_blkdiagICBC(obj,x0,k,kind) % for Tfull
        x = apply_scpc_fourier_Vav_blkdiagICBC(obj,x0,k,kind) % for T
        x = apply_Tfull_star(obj,x0,k);
        x = apply_Tfull_star_Tfull(obj,x0,k);
        
        % default preconditioner for T
        function x = apply_pc(obj,x0,k)
            x = obj.apply_scpc_fourier_Vav_blkdiagICBC(x0,k,'inv');
        end
        
        % default preconditioner for Tfull
        function x = apply_pc_full(obj,x0,k)
            x = obj.apply_pc_fourier_Vav_blkdiagICBC(x0,k,'inv');
        end
        
        function x = apply_fourierA12(obj,x0)
           % x = fourierA12*x0
           wx0 = zeros(length(obj.ratf),1);
           for j = 1:length(obj.ratf)
               xj = x0(obj.ratf_ends(j)-obj.ratf_lens(j)+1 : obj.ratf_ends(j));
               wx0(j) = obj.ratf(j).w*xj;
           end
           x = zeros(obj.dtns.Nr*obj.dtns.Nt,1); x(obj.dtns.BCrows) = wx0;
        end
            
        function x = apply_fourierA12_star(obj,x0)
            % x = fourierA12'*x0
	    x = zeros(sum(obj.ratf_lens),1);
            for j = 1:length(obj.ratf_lens)
		x(obj.ratf_ends(j)-obj.ratf_lens(j)+1 : obj.ratf_ends(j)) = ...
                    obj.ratf(j).w'*x0(obj.dtns.BCrows(j));
            end
        end

        function x = apply_A12_star(obj,x0)
	    % x = A12'*x0
            x = obj.dtns.Nt*obj.apply_fourierA12_star(obj.dtns.apply_Uinv_kron_I(x0));
        end

        function x = apply_fourierA21(obj,x0)
            % x = fourierA21*x0
            x = zeros(sum(obj.ratf_lens),1);
            for j = 1:length(obj.ratf_lens)
                x(obj.ratf_ends(j)-obj.ratf_lens(j)+1 : obj.ratf_ends(j)) = ...
                    x0(obj.dtns.BCrows(j));
            end
        end

        function x = apply_fourierA21_star(obj,x0)
            % x = fourierA21'*x0
           onesx0 = zeros(length(obj.ratf),1);
           for j = 1:length(obj.ratf)
               xj = x0(obj.ratf_ends(j)-obj.ratf_lens(j)+1 : obj.ratf_ends(j));
               onesx0(j) = sum(xj);
           end
           x = zeros(obj.dtns.Nr*obj.dtns.Nt,1); x(obj.dtns.BCrows) = onesx0;
        end

        function x = apply_A21_star(obj,x0)
	    % x = A21'*x0
            x = obj.dtns.apply_U_kron_I(obj.apply_fourierA21_star(x0))/obj.dtns.Nt;
        end
    end
    methods (Static)
        ratApproxPoles = get_poles_in_range(n,region,R,br)
    end
end

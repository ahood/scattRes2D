classdef DtNBC_sparse < scattResComp2d_sparse
    % given a smooth potential (at least smooth on each sub-annulus) 
    % with compact support, discretizes and solves
    % scattering problem with DtN boundary conditions.    
    properties
        VvaluesVec % vector of potential values along rays
        coords
        A_blks_common % part of blocks of A that is same for all blocks
        Vhat % VvaluesVec but reshaped
        Vav % average of potential on concentric circles
        blk_precond_k % the value of k last used for preconditioner
        M_blks % a cell array where j-th element is struct with L,U decomp of j-th block of block diag preconditioner M  
        blk_fourier_precond_k
        fourierM_blks % a cell array j-th elt is struct with L,U decomp of j-th block
    end
    methods
        function obj = DtNBC_sparse(Nt,Nrs,Vs,coords,Rs)
        % V is a potential function handle. V must be smooth on each
        % submesh.
            obj@scattResComp2d_sparse(Nt,Nrs,Rs);
            if isnumeric(Vs)
                obj.VvaluesVec = Vs;
            else
                obj.VvaluesVec = obj.valuesVecFromFunCellArray(Vs,coords);
            end
            obj.coords = coords;
            obj.ICBCmatrices();
            
            % for apply_pc_blkdiag
            obj.A_blks_common = -diag(1./obj.r)*obj.Dr_I - obj.Drr_I - obj.Dtt(1,1)*diag(1./obj.r.^2);
            % for apply_pc_blkdiag, apply_pc_diag and fourier_V_blkdiagICBC
            obj.Vhat = reshape(obj.VvaluesVec, obj.Nr, []);
            
            % for fourier_Vav_blkdiagICBC
            obj.Vav = mean(obj.Vhat,2);
            
            % for apply_pc_blkdiag
            obj.blk_precond_k = Inf;
            obj.M_blks = cell(1,obj.Nt);

            % for apply_pc_fourier_Vav_blkdiagICBC -- current default
            obj.blk_fourier_precond_k = Inf;
            obj.fourierM_blks = cell(1,obj.Nt);
        end
        
        % don't form, apply
        x = apply_T(obj,x0,k) % x = T(k)*x0
        x = apply_pc_blkdiag(obj,x0,k) % x = M(k)\x0, M(k) the blockdiag of T(k)
        x = apply_pc_fourier_Vav_blkdiagICBC(obj,x0,k,kind) % use block diag part of ICBCs in Fourier space
        x = apply_A_star(obj,x0); % x = A'*x0
        
        % default preconditioner
        function x = apply_pc(obj,x0,k)
            x = obj.apply_pc_fourier_Vav_blkdiagICBC(x0,k,'inv');
        end
    end
    methods (Static)
        % DtN map and its derivatives w.r.t. k
        f = DtNcoeffs(n,k,R)
        df = dDtNcoeffs(n,k,R)
        d2f = d2DtNcoeffs(n,k,R)
    end
end

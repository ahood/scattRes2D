classdef dirBC_sparse < scattResComp2d_sparse
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with Dirichlet boundary conditions.
    
    properties
        VvaluesVec
        coords
    end
    methods
        function obj = dirBC_sparse(Nt,Nrs,Vs,coords,Rs)
        % V is a potential function handle. V must be smooth on the R
        % disk.
            obj@scattResComp2d_sparse(Nt,Nrs,Rs);
            if isnumeric(Vs)
                obj.VvaluesVec = Vs;
            else
                obj.VvaluesVec = obj.valuesVecFromFunCellArray(Vs,coords);
            end
            obj.coords = coords;
        end
        
        % compute E such that A - E*B singular
        [evecs,evals] = eig_comp(obj,n,E0,verbose)
        
        % default preconditioner - identity for now
        function x = apply_pc(obj,x0,k)
            x = x0; 
        end

        % don't form, apply
        x = apply_T(obj,x0,k)
        x = apply_B(obj,x0)
    end
end
classdef dirBC < scattResComp2d
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with Dirichlet boundary conditions.
    
    properties
        A      
        B      % A-E*B singular characterizes approximate resonances as sqrt(E) values
        Aorig
        Borig % A and B before boundary and interface conditions applied
        A_fc
        B_fc
        A_cc
        B_cc % A and B in Chebyshev basis
        Vs
        coords % coordinates of the potential function (rect,polar,complex)
        Lr
        Lr_unscaled
        evals
        evecs
    end
    methods
        function obj = dirBC(Nt,Nrs,Vs,coords,Rs)
        % V is a potential function handle. V must be smooth on the R
        % disk.
            obj@scattResComp2d(Nt,Nrs,Rs);
            
            obj.Vs = Vs;
            obj.coords = coords;
            obj.Lr  = diag(1./obj.rs)*obj.Dr  + obj.Drr; % radial derivatives part of Laplacian
            obj.enforceDE(obj.Lr,Vs,coords,obj.r,obj.r);
            obj.enforceIC();
            obj.enforceBC();
            
            % fourier basis
            obj.enforceDE_fc();
            obj.enforceIC_fc();
            obj.enforceBC_fc();

            % cheb basis
            obj.enforceDE_cc();
            obj.enforceIC_cc();
            obj.enforceBC_cc();
        end
        function eig_comp(obj)
            [evecs,evals] = eig(full(obj.A),full(obj.B));
            [evecs,evals] = obj.sort_eigs(evecs,evals);
            obj.evals = evals; obj.evecs = evecs;            
        end
        function enforceBC(obj)
            obj.A(obj.BCrows,:) = 0; obj.B(obj.BCrows,:) = 0;
            for row = obj.BCrows, obj.A(row,row) = 1; end            
        end        
        function Tk = T(obj,k)
            Tk = obj.A - k^2*obj.B;
        end
        function Tk = T_fc(obj,k)
            Tk = obj.A_fc - k^2*obj.B_fc;
        end
        function Tk = T_cc(obj,k)
            Tk = obj.A_cc - k^2*obj.B_cc;
        end
        
        % enforce boundary conditions
        enforceBC_fc(obj)
        enforceBC_cc(obj)
    end
end
classdef dirBC_axisymm < scattResComp2d_axisymm
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with Dirichlet boundary conditions.
    
    properties
        A_fc_blks
        B_fc_blks % cell arrays of diagonal blocks in Fourier basis
        T_fc_blks
        A_cc_blks
        B_cc_blks % cell arrays of diagonal blocks in Cheb basis
        T_cc_blks
        Vs
        coords % coordinates of the potential function (rect,polar,complex)
    end
    methods
        function obj = dirBC_axisymm(Nt,Nrs,Vs,coords,Rs)
        % V is a potential function handle. V must be smooth on the R
        % disk.
            obj@scattResComp2d_axisymm(Nt,Nrs,Rs);
            
            obj.Vs = Vs;
            obj.coords = coords;

            % fourier basis
            obj.enforceDE_fc();
            obj.enforceIC_fc();
            obj.enforceBC_fc();

            % cheb basis
            obj.enforceDE_cc();
            obj.enforceIC_cc();
            obj.enforceBC_cc();
            
            obj.Es = cell(1,Nt);
            obj.ks = cell(1,Nt);
        end
        function T_fc(obj,k)
            for j = 1:obj.Nt
                obj.T_fc_blks{j} = obj.A_fc_blks{j} - k^2*obj.B_fc_blks{j};
            end
        end
        function T_cc(obj,k)
            for j = 1:obj.Nt
                obj.T_cc_blks{j} = obj.A_cc_blks{j} - k^2*obj.B_cc_blks{j};
            end
        end
        function Tk = T_fc_j(obj,k,j)
            Tk = obj.A_fc_blks{j} - k^2*obj.B_fc_blks{j};
        end
        function Tk = T_cc_j(obj,k,j)
            Tk = obj.A_cc_blks{j} - k^2*obj.B_cc_blks{j};
        end
        
        % enforce boundary conditions
        enforceBC_fc(obj)
        enforceBC_cc(obj)
    end
end
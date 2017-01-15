classdef DtNBC_axisymm < scattResComp2d_axisymm
    % Given a smooth potential (at least smooth on each sub-annulus) 
    % with compact support, discretizes and solves
    % scattering problem with DtN boundary conditions.    
    properties
        A_fc_blks
        B_fc_blks  
        C_fc_blks  
        T_fc_blks  % cell array of diagonal blocks in Fourier basis
        A_cc_blks
        B_cc_blks  
        C_cc_blks  
        T_cc_blks  % same, except in Cheb coeff basis
        dC_fc_blks
        Vs     % cell array of potential functions for each annular region
        coords % coordinates of the potential function (rect,polar,complex)
    end
    methods
        function obj = DtNBC_axisymm(Nt,Nrs,Vs,coords,Rs)
        % V is a potential function handle. V must be smooth on each
        % submesh.
            obj@scattResComp2d_axisymm(Nt,Nrs,Rs);

            obj.Vs = Vs;
            obj.coords = coords;

            % Fourier basis
            obj.enforceDE_fc();
            obj.enforceIC_fc();
            obj.enforceBC_fc();

            % Cheb basis
            obj.enforceDE_cc();
            obj.enforceIC_cc();
            obj.enforceBC_cc();
        end
        
        % enforce boundary conditions
        enforceBC_fc(obj) % Fourier basis
        enforceBC_cc(obj) % Cheb coeff basis
        
        function T_fc(obj,k)
            obj.C_fc(k);
            for j = 1:obj.Nt
                obj.T_fc_blks{j} = obj.A_fc_blks{j} - k^2*obj.B_fc_blks{j} + ...
                    obj.C_fc_blks{j};
            end
        end
        function T_cc(obj,k)
            obj.C_cc(k);
            for j = 1:obj.Nt
                obj.T_cc_blks{j} = obj.A_cc_blks{j} - k^2*obj.B_cc_blks{j} + ...
                    obj.C_cc_blks{j};
            end
        end
        function Ck = C_fc_j(obj,k,j)
            n = obj.Ns(j);
            Ck = spalloc(obj.Nr,obj.Nr,1);
            Ck(end,end) = obj.DtNcoeffs(n,k,obj.r(end));
        end
        function dC_fc(obj,k)
            d = obj.dDtNcoeffs(obj.Ns,k,obj.r(end));
            for j = 1:obj.Nt
                obj.dC_fc_blks{j} = spalloc(obj.Nr,obj.Nr,1);
                obj.dC_fc_blks{j}(end,end) = d(j);
            end
        end
        function Tk = T_fc_j(obj,k,j)
            Tk = obj.A_fc_blks{j} - k^2*obj.B_fc_blks{j};
            n = obj.Ns(j);
            Tk(end,end) = Tk(end,end) + obj.DtNcoeffs(n,k,obj.r(end));
        end
        function dTk = dT_fc_j(obj,k,j)
            dTk = -2*k*obj.B_fc_blks{j};
            n = obj.Ns(j);
            dTk(end,end) = dTk(end,end) + obj.dDtNcoeffs(n,k,obj.r(end));
        end
        
        C_fc(obj,k)
        C_cc(obj,k)
    end
    methods (Static)
        % DtN map and its derivatives w.r.t. k
        f = DtNcoeffs(n,k,R)
        df = dDtNcoeffs(n,k,R)
        d2f = d2DtNcoeffs(n,k,R)
    end
end
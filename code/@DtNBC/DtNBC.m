classdef DtNBC < scattResComp2d
    % Given a smooth potential (at least smooth on each sub-annulus) 
    % with compact support, discretizes and solves
    % scattering problem with DtN boundary conditions.    
    properties
        Aorig % A before BCs applied, too
        Borig
        A    
        B     % eigs of A + (DtN thing depending on k) - k^2*B characterizes resonances
        A_fc
        B_fc  % same, except in Fourier basis
        A_cc
        B_cc  % same, except in Cheb coeff basis
        Vs    % cell array of potential functions for each annular region
        coords % coordinates of the potential function (rect,polar,complex)
        Lr
    end
    methods
        function obj = DtNBC(Nt,Nrs,Vs,coords,Rs)
        % V is a potential function handle. V must be smooth on each
        % submesh.
            obj@scattResComp2d(Nt,Nrs,Rs);

            obj.Vs = Vs;
            obj.coords = coords;
            obj.Lr = sparse(diag(1./obj.rs))*obj.Dr  + obj.Drr; % radial derivatives part of Laplacian
            obj.enforceDE(obj.Lr,Vs,coords,obj.r,obj.r);

            % grab rows that are about to be replaced
            obj.ArowChange = -obj.A(obj.ICBCrows,:);
            obj.BrowChange = -obj.B(obj.ICBCrows,:);

            % change A and B
            obj.enforceIC();
            obj.enforceBC();

            % Fourier basis
            obj.enforceDE_fc();
            obj.enforceIC_fc();
            obj.enforceBC_fc();

            % Cheb basis
            obj.enforceDE_cc();
            obj.enforceIC_cc();
            obj.enforceBC_cc();

            % incorporate new rows
            obj.ArowChange = obj.ArowChange + obj.A(obj.ICBCrows,:);
            obj.BrowChange = obj.BrowChange + obj.B(obj.ICBCrows,:);
        end
        
        % enforce boundary conditions
        enforceBC(obj)
        enforceBC_fc(obj) % Fourier basis
        enforceBC_cc(obj) % Cheb coeff basis
        
        % computes matrix-valued function T whose eigs are the resonances
        Tk = T(obj,k)
        dTk = dT(obj,k) % T'(k)
        d2Tk = d2T(obj,k) % T''(k)

        % Fourier basis
        Tk = T_fc(obj,k)
        
        % Cheb coeff basis
        Tk = T_cc(obj,k) 
        
        % nonpolynomial part C(k) of T(k)
        Ck = C(obj,k)
        Ck = C_fc(obj,k)
        Ck = C_cc(obj,k)
        
        % derivative of Log(T(sqrt(z))), branch given by br, for winding #
        dfoverf = dLogTz(obj,z,br)
        d_dfoverf = d2LogTz(obj,z,br)
    end
    methods (Static)
        % DtN map and its derivatives w.r.t. k
        f = DtNcoeffs(n,k,R)
        df = dDtNcoeffs(n,k,R)
        d2f = d2DtNcoeffs(n,k,R)
    end
end
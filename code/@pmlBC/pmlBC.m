classdef pmlBC < scattResComp2d
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with a PML and Dirichlet boundary 
    % conditions.
    
    properties
        s
        ss
        Ds
        Dss
        Aorig
        Borig
        A      
        B      % A-E*B singular characterizes approximate resonances as sqrt(E) values
        A_fc
        B_fc   % Fourier basis
        A_cc
        B_cc   % Cheb basis
        Vs     % cell array of potentials on annular regions
        coords % coordinates of the potential function (rect,polar,complex)
        Lr
        pieces % schurComp object storing indices which do and do not correspond to DtN problem
        evals
        evecs
        ds
        d2s
        ds_short
        d2s_short
    end
    methods
        function obj = pmlBC(Nt,Nrs_in,Nr_out,Vs,coords,Rs_in,R_out,l,dl,d2l)
        % Rout is the place Dirichlet BCs will be applied. The PML variable
        % satisfies s(r) = r + 1i*l(r) for r > Rin, and dl
        % and d2l are derivatives of l.
            Nrs = Nrs_in; Nrs(end + 1) = Nr_out;
            Rs  = Rs_in;  Rs (end + 1) = R_out;
            obj@scattResComp2d(Nt,Nrs,Rs);
            
            obj.Vs = Vs;
            obj.coords = coords;
            
            R_in = Rs_in(end);
            r2s   = @(t) t + (t > R_in).*1i.*  l(t);
            dr2s  = @(t) 1 + (t > R_in).*1i.* dl(t);
            d2r2s = @(t)     (t > R_in).*1i.*d2l(t);

            obj.s = r2s(obj.r);
            obj.ss = r2s(obj.rs);
            obj.ds = dr2s(obj.rs); 
            obj.d2s = d2r2s(obj.rs);
            obj.ds_short = dr2s(obj.r);
            obj.d2s_short = d2r2s(obj.r);
            
            obj.Ds = sparse(diag(1./obj.ds))*obj.Dr;
            obj.Dss = -sparse(diag(obj.d2s./obj.ds.^3))*obj.Dr + ...
                       sparse(diag(1./obj.ds.^2))*obj.Drr;   

            obj.Lr  = sparse(diag(1./obj.ss))*obj.Ds  + obj.Dss; % radial derivatives part of Laplacian
            obj.enforceDE(obj.Lr,Vs,coords,obj.r,obj.s);
            obj.enforceIC();
            obj.enforceBC();
 
            Nr_in = sum(Nrs_in);
            I1 = kron(ones(obj.Nt,1),(1:Nr_in)') + kron((0:obj.Nt-1)',0*(1:Nr_in)' + obj.Nr);
            I2 = 1:obj.Nr*obj.Nt; I2(I1) = [];
            obj.pieces = schurComp(I1,I2);

            % Fourier basis
            obj.enforceDE_fc();
            obj.enforceIC_fc();
            obj.enforceBC_fc();

            % Cheb basis
            obj.enforceDE_cc();
            obj.enforceIC_cc();
            obj.enforceBC_cc();
        end
        function eig_comp(obj)
            % compute eigenpairs and sort by magnitude
            [evecs,evals] = eig(full(obj.A),full(obj.B));
            [evecs,evals] = obj.sort_eigs(evecs,evals);
            obj.evals = evals; obj.evecs = evecs;            
        end        
        function enforceBC(obj)    
            % dirichlet boundary condition
            obj.A(obj.BCrows,:) = 0; obj.B(obj.BCrows,:) = 0;
            for row = obj.BCrows, obj.A(row,row) = 1; end
        end
        
        % solve the scattering problem on disc including PML region
        scattValuesVec_ext = solve_full(obj,k,incfun)
        scattValuesVec_ext = solve_full_fc(obj,k,incfun)
        scattValuesVec_ext = solve_full_cc(obj,k,incfun)

        % RHS for full and schur complemented scattering problems
        RHS = RHSfromFun(obj,incfun)
        RHS = RHSfromFun_fc(obj,incfun)
        RHS = RHSfromFun_cc(obj,incfun)
        RHS = RHSfromFun_full(obj,incfun)
        RHS = RHSfromFun_full_fc(obj,incfun)
        RHS = RHSfromFun_full_cc(obj,incfun)
        
        % schur complement which approximates DtNBC.T(k)
        Tk  = T(obj,k)
        dTk = dT(obj,k)
        Tk = Tfull(obj,k) % Tk = A - k^2*B

        % in Fourier and Cheb bases
        function Tk = Tfull_fc(obj,k)
            Tk = obj.A_fc - k^2*obj.B_fc; 
        end
        function Tk = Tfull_cc(obj,k)
            Tk = obj.A_cc - k^2*obj.B_cc;
        end
        function Tk = T_fc(obj,k)
            Tk = obj.pieces.S(obj.A_fc - k^2*obj.B_fc);
        end
        function Tk = T_cc(obj,k)
            Tk = obj.pieces.S(obj.A_cc - k^2*obj.B_cc);
        end
        
        % its nonpolynomial part
        Ck = C(obj,k)

        % enforcing BCs in Fourier and Cheb bases
        enforceBC_fc(obj)
        enforceBC_cc(obj)
    end
end    
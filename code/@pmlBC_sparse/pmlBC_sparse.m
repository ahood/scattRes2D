classdef pmlBC_sparse < scattResComp2d_sparse
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with a PML and Dirichlet boundary 
    % conditions.
    
    properties
        s
        ss
        ds
        d2s
        VvaluesVec
        coords
        pieces % schurComp object storing indices which do and do not correspond to DtN problem
    end
    methods
        function obj = pmlBC_sparse(Nt,Nrs_in,Nr_out,Vs,coords,Rs_in,R_out,l,dl,d2l)
        % Rout is the place Dirichlet BCs will be applied. The PML variable
        % satisfies s(r) = r + 1i*l(r) for r > Rin, and dl
        % and d2l are derivatives of l.
            Nrs = Nrs_in; Nrs(end + 1) = Nr_out;
            Rs  = Rs_in;  Rs (end + 1) = R_out;
            obj@scattResComp2d_sparse(Nt,Nrs,Rs);
            if isnumeric(Vs)
                obj.VvaluesVec = Vs;
            else
                obj.VvaluesVec = obj.valuesVecFromFunCellArray(Vs,coords);
            end
            obj.coords = coords;
            
            R_in = Rs_in(end);
            r2s   = @(t) t + (t > R_in).*1i.*  l(t);
            dr2s  = @(t) 1 + (t > R_in).*1i.* dl(t);
            d2r2s = @(t)     (t > R_in).*1i.*d2l(t);

            obj.s = r2s(obj.r);
            obj.ss = r2s(obj.rs);
            obj.ds = dr2s(obj.rs); 
            obj.d2s = d2r2s(obj.rs);
            
            obj.R = spdiags(1./obj.s,0,obj.Nr,obj.Nr);
 
            Nr_in = sum(Nrs_in);
            I1 = kron(ones(obj.Nt,1),(1:Nr_in)') + kron((0:obj.Nt-1)',0*(1:Nr_in)' + obj.Nr);
            I2 = 1:obj.Nr*obj.Nt; I2(I1) = [];
            obj.pieces = schurComp(I1,I2);
        end
        
        % compute resonances iteratively
        [evecs,evals] = eig_comp(obj,n,E0,verbose)
        
        % default preconditioner - identity for now
        function x = apply_pc_full(obj,x0,k)
            x = x0; 
        end
        
        % solve the scattering problem on disc including PML region
        scattValuesVec_ext = solve_full(obj,k,incfun,precond,verbose)
        
        % don't form, apply
        x = apply_op_no_bcs(obj,x0,E0,VvaluesVec) % override
        x = apply_Tfull(obj,x0,k) % x = Tfull(k)*x0
        x = apply_Bfull(obj,x0)
        
        % RHS for scattering problems
        RHS = RHSfromFun(obj,incfun)
        RHS = RHSfromFun_full(obj,incfun)
    end
end    
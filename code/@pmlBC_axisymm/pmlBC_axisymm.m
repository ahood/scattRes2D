classdef pmlBC_axisymm < scattResComp2d_axisymm
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with a PML and Dirichlet boundary 
    % conditions.
    
    properties
        Nr_in
        s
        ss
        Afull_fc_blks
        Bfull_fc_blks   
        Tfull_fc_blks   % Fourier basis
        Afull_cc_blks
        Bfull_cc_blks
        Tfull_cc_blks   % Cheb basis
        A_fc_blks
        B_fc_blks
        C_fc_blks
        T_fc_blks
        A_cc_blks
        B_cc_blks
        C_cc_blks
        T_cc_blks
        Vs     % cell array of potentials on annular regions
        coords % coordinates of the potential function (rect,polar,complex)
        pieces % schurComp object storing indices which do and do not correspond to DtN problem
        ds
        d2s
        ds_short
        d2s_short
    end
    methods
        function obj = pmlBC_axisymm(Nt,Nrs_in,Nr_out,Vs,coords,Rs_in,R_out,l,dl,d2l)
        % Rout is the place Dirichlet BCs will be applied. The PML variable
        % satisfies s(r) = r + 1i*l(r) for r > Rin, and dl
        % and d2l are derivatives of l.
            Nrs = Nrs_in; Nrs(end + 1) = Nr_out;
            Rs  = Rs_in;  Rs (end + 1) = R_out;
            obj@scattResComp2d_axisymm(Nt,Nrs,Rs);
            
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
            
            Nr_in = sum(Nrs_in);
            obj.Nr_in;
            I1 = 1:Nr_in;
            I2 = 1:obj.Nr; I2(I1) = [];
            obj.pieces = schurComp(I1,I2);

            % Fourier basis
            obj.enforceDE_full_fc();
            obj.enforceIC_full_fc();
            obj.enforceBC_full_fc();
            for j = 1:obj.Nt
                obj.A_fc_blks{j} = obj.Afull_fc_blks{j}(1:Nr_in,1:Nr_in);
                obj.B_fc_blks{j} = obj.Bfull_fc_blks{j}(1:Nr_in,1:Nr_in);
            end
            % will want to put derivative then continuity condition
            % for last interface if we want to interpret schur complement
            % part as DtN map approx

            % Cheb basis
            obj.enforceDE_full_cc();
            obj.enforceIC_full_cc();
            obj.enforceBC_full_cc();
            for j = 1:obj.Nt
                obj.A_cc_blks{j} = obj.Afull_cc_blks{j}(1:Nr_in,1:Nr_in);
                obj.B_cc_blks{j} = obj.Bfull_cc_blks{j}(1:Nr_in,1:Nr_in);
            end
        end
        
        % solve the scattering problem on disc including PML region
        scattValuesVec_ext = solve_full_fc(obj,k,incfun)
        scattValuesVec_ext = solve_full_cc(obj,k,incfun)

        % RHS for full and schur complemented scattering problems
        RHS = RHSfromFun_fc(obj,incfun)
        RHS = RHSfromFun_cc(obj,incfun)
        RHS = RHSfromFun_full_fc(obj,incfun)
        RHS = RHSfromFun_full_cc(obj,incfun)
        
        % T(k) in Fourier and Cheb bases
        function Tfull_fc(obj,k)
            for j = 1:obj.Nt
                obj.Tfull_fc_blks{j} = obj.Afull_fc_blks{j} - k^2*obj.Bfull_fc_blks{j};
            end
        end
        function Tfull_cc(obj,k)
            for j = 1:obj.Nt
                obj.Tfull_cc_blks{j} = obj.Afull_cc_blks{j} - k^2*obj.Bfull_cc_blks{j};
            end
        end
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
        
        % enforcing BCs in Fourier and Cheb bases
        C_fc(obj,k)
        C_cc(obj,k)
        enforceDE_full_fc(obj)
        enforceIC_full_fc(obj)
        enforceBC_full_fc(obj)
        enforceDE_full_cc(obj)
        enforceIC_full_cc(obj)
        enforceBC_full_cc(obj)
    end
end    
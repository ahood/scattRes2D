classdef ratApproxDtNBC_axisymm < handle
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with approximate DtN boundary conditions,
    % where approximations are meant to be good in region
    properties
        Afull_fc_blks
        Bfull_fc_blks
        Tfull_fc_blks
        C_fc_blks
        T_fc_blks 
        Aschur
        
        mysqrt
        Es
        ks
        br = 1; % branch of sqrt (1 or 2)
        dtn % DtNBC_axisymm problem (I think this will just be a reference since subclass of handle, so no memory issue?)
        ell % ellipse in which the rational approximation is supposed to be pretty good
        ratf = []; % list of rational approx objects for individual components of DtN map
        poles = {};
        pieces
    end
    methods
        function obj = ratApproxDtNBC_axisymm(dtn,ell,br,N)
        % V is a potential function handle. V must be smooth on the R
        % disk and V = V(x,y), not V = V(r,theta).
            if nargin < 4, N = 500; end
            obj.br = br; % 1 or 2
            obj.dtn = dtn;
            obj.ell = ell;
            
            % use given branch of square root function
            if br == 1, mysqrt = @(z) sqrt(z);
            elseif br == 2, mysqrt = @(z) -sqrt(z);
            end
            obj.mysqrt = mysqrt;
            
            for j = 1:dtn.Nt
                n = dtn.Ns(j);
                fn = @(z) DtNBC.DtNcoeffs(n,mysqrt(z),dtn.r(end));
                poles{j} = ratApproxDtNBC.get_poles_in_range(n,ell,dtn.r(end),br);
                t = linspace(0,1,N);
                rs = 0*poles{j} + 0.05; % wow, also so dumb
                ratfn = ratApproxWithPoles(ell,fn,poles{j},rs,t);
                ratf(j) = ratfn;
                
                % make Afull and Bfull blocks
                obj.Afull_fc_blks{j} = blkdiag(dtn.A_fc_blks{j}, diag(ratfn.z));
                obj.Afull_fc_blks{j}(dtn.Nr,dtn.Nr+1:end) = ratfn.w;
%                 obj.Afull_fc_blks{j}(dtn.Nr+1:end,dtn.Nr) = 0*ratfn.w' + 1;
                obj.Afull_fc_blks{j}(dtn.Nr+1:end,dtn.Nr) = 0*ratfn.w' - 1;
                obj.Bfull_fc_blks{j} = blkdiag(dtn.B_fc_blks{j}, eye(length(ratfn.z)));
            end
            obj.ratf = ratf;
            obj.poles = poles;
                
            % make Aschur
            I2 = sort([dtn.derRows(:); dtn.valRows(:); dtn.Nr]); % rows where Bfull is zero (same for all)
            for j = 1:dtn.Nt
                I1 = 1:length(obj.Afull_fc_blks{j}); I1(I2) = [];
                sc = schurComp(I1,I2);
                obj.Aschur{j} = sc.S(obj.Afull_fc_blks{j});
            end
            
            I1 = 1:dtn.Nr; 
            for j = 1:dtn.Nt
                n = length(ratf(j).z);
                I2 = dtn.Nr + (1:n);
                obj.pieces{j} = schurComp(I1,I2);
            end
        end
        function Tfull_fc(obj,k)
            for j = 1:obj.dtn.Nt
                obj.Tfull_fc_blks{j} = obj.Afull_fc_blks{j} - k^2*obj.Bfull_fc_blks{j};
            end
        end
        function Tk = Tfull_fc_j(obj,k,j)
            Tk = obj.Afull_fc_blks{j} - k^2*obj.Bfull_fc_blks{j};
        end
        function T_fc(obj,k)
            obj.C_fc(k);
            for j = 1:obj.dtn.Nt
                obj.T_fc_blks{j} = obj.dtn.A_fc_blks{j} ...
                             - k^2*obj.dtn.B_fc_blks{j} ...
                             + obj.C_fc_blks{j};
            end
        end
        function Tk = T_fc_j(obj,k,j)
            Tk = obj.dtn.A_fc_blks{j} - k^2*obj.dtn.B_fc_blks{j};
            Tk(end,end) = Tk(end,end) + obj.ratf(j).eval_at(k^2);
        end
        function Ck = C_fc_j(obj,k,j)
            Ck = spalloc(obj.dtn.Nr, obj.dtn.Nr,1);
            Ck(end,end) = obj.ratf(j).eval_at(k^2);
        end
        function C_fc(obj,k)
            z = k^2;
            for j = 1:obj.dtn.Nt
                obj.C_fc_blks{j} = spalloc(obj.dtn.Nr,obj.dtn.Nr,1);
                obj.C_fc_blks{j}(end,end) = obj.ratf(j).eval_at(z);
            end
        end
        function eig_comp(obj)
            for j = 1:obj.dtn.Nt
%                 obj.Es{j} = eig(obj.Afull_fc_blks{j},full(obj.Bfull_fc_blks{j}));
%                 obj.Es{j} = obj.Es{j}( ~isinf(obj.Es{j}) );
                obj.Es{j} = eig(obj.Aschur{j});
                obj.ks{j} = obj.mysqrt(obj.Es{j});
            end
        end

        % solve the scattering problem
        scattValuesVec = solve(obj,k,incfun)
        
        % look at error in rational approximation on given rectangle r
        E = compute_error(obj,r)
        show_error(obj,r)
        
        % look at poles
        show_poles(obj,r)
        
        % for scattering with incident wave exp(i k x)
        Tk = T(obj,k)     % schur complement of A - k^2*B
        Tk = Tfull(obj,k) % plain old A - k^2*B
    end
    methods (Static)
        ratApproxPoles = get_poles_in_range(n,region,R,br)
    end
end
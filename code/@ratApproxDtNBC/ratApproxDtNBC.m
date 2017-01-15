classdef ratApproxDtNBC < handle
    % given a smooth potential with compact support, discretizes and solves
    % scattering and resonance problems with approximate DtN boundary conditions,
    % where approximations are meant to be good in region
    properties
        A   
        B % matrices whose eigenvalues approximate resonance energies (E = k^2)
        A12
        A21
        fourierA12
        fourierA21
        Aschur % matrix with same eigs as (A,B) pencil
        br = 1; % branch of sqrt (1 or 2)
        dtn % DtN problem (I think this will just be a reference since subclass of handle, so no memory issue?)
        ell % ellipse in which the rational approximation is supposed to be pretty good
        ratf = []; % list of rational approx objects for individual components of DtN map
        poles = {};
        sc % schur complement object for map from scattered wave to -V*incident
        evals
        evecs
    end
    methods
        function obj = ratApproxDtNBC(dtn,ell,br,N)
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
%             ns = -dtn.Ns:dtn.Ns;
            Idtn = (1:dtn.Nt)*dtn.Nr; % indices where DtN map components go
            % set up rational approximation to dtn map part of
            % dtn.T(mysqrt(z))
            A11 = dtn.A; fourierA12 = []; fourierA21 = []; A22diag = [];
            B11 = dtn.B;
            
%             for jj = length(ns):-1:1 % doing twice the work here
%             for jj = 1:length(ns) % I want to preallocate array of objects but I don't know how to do so nicely if the object needs args
%                 n = ns(jj);
            for jj = 1:dtn.Nt
                n = dtn.Ns(jj);
                fn = @(z) DtNBC.DtNcoeffs(n,mysqrt(z),dtn.r(end));
                poles{jj} = ratApproxDtNBC.get_poles_in_range(n,ell,dtn.r(end),br);
                t = linspace(0,1,N);
                rs = 0*poles{jj} + 0.05; % wow, also so dumb TODO: make the user choose
                ratfn = ratApproxWithPoles(ell,fn,poles{jj},rs,t);
                ratf(jj) = ratfn;
                
                % append onto A and B
                row = Idtn(jj);
                cols = size(fourierA12,2) + (1:length(ratfn.w));
                fourierA12(row ,cols) = ratfn.w;
                fourierA21(cols,row ) = 0*ratfn.w' + 1;
                A22diag( cols) = ratfn.z;
            end
            obj.ratf = ratf;
            obj.poles = poles;
            
            % Account for change of basis
            A12 =            kron(obj.dtn.U   , eye(obj.dtn.Nr))*fourierA12;
            A21 = fourierA21*kron(obj.dtn.Uinv, eye(obj.dtn.Nr));
            obj.A12 = A12; obj.A21 = A21;

            obj.fourierA12 = fourierA12; obj.fourierA21 = fourierA21;
            obj.A = sparse([A11, A12; A21, diag(A22diag)]);
            obj.B = sparse(blkdiag(B11, eye(length(A22diag)) ));
            
            I1 = find( diag(obj.B) ); % get indices where 1s on diag
            I2 = 1:length(obj.B); I2(I1) = [];
            ratschur = schurComp(I1,I2);
            obj.Aschur = ratschur.S(obj.A);
            obj.sc = schurComp(1:length(A11), length(A11)+(1:length(A22diag)));
        end
        
        % compute n resonances closest to E0
        function [resvecs,resvals] = resonances(obj,n,E0)
            % Uses eigs to find the n eigs of (A,B) closest to E0
            if nargin < 3, E0 = obj.ell.c; end
            if nargin < 2 || isempty(n), n = 40; end
            opts.isreal = false;
            [resvecs,resvals] = eigs(obj.A, obj.B, n, E0, opts);
            resvals = diag(resvals);
        end
        function eig_comp(obj,region)
            [evecs,evals] = eig(full(obj.A),full(obj.B));
            if nargin == 2
                idx = region.contains(evals);
                evals = evals(idx); evecs = evecs(:,idx);
            end
            [evecs,evals] = obj.dtn.sort_eigs(evecs,evals);
            obj.evals = evals; obj.evecs = evecs;            
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
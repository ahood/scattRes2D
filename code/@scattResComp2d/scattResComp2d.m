classdef scattResComp2d < handle
    % General framework for scattering computation and plotting.
    properties
        Nt     % number of thetas (must be even)
        Ns
        Nrs    % size of meshes between breakpoints
        Nr     % total number of mesh points (including duplicates)
        Rs     % breakpoints for r direction
        r      % mesh in r direction
        Dr     % mapping from values to derivatives in r at mesh points
        Drr    % ditto, but second derivatives in r
        Dtt    % second derivatives in theta
        Drs_even % maps radial values of even fcn to radial values of its deriv (which is odd )
        Drs_odd  % maps radial values of odd  fcn to radial values of its deriv (which is even)
        Drs_even_unscaled
        Drs_odd_unscaled
        theta  % mesh of angles
        rr
        tt     % 2D mesh in r,theta coords
        xx
        yy     % 2D mesh in x,y coords
        xs
        ys     % x,y coords used for valuesVec
        rs
        ts     % r,theta coordinates used for valuesVec
        zs     % complex coord used for values Vec
        U
        Uinv   % needed for applying DtN map as pseudodifferential op
        Winv_even
        Winv_odd
        Dr_cc
        Drs_cc_even % maps cheb coeffs of even function to cheb coeffs of its deriv (an odd  fcn)
        Drs_cc_odd  % maps cheb coeffs of odd  function to cheb coeffs of its deriv (an even fcn)
        BCrows     % rows corresponding to boundary of 2d region
        valRows    % rows in each subblock for continuity at interface (Fourier basis)
        derRows    %                           C1
        valRows_cc %         (Cheb basis)
        derRows_cc %         (Cheb basis)
        ICBCrows % the locations of all rows changed for ICs or BCs
        
        % stuff I need to make fast matvec
        Dr_I
        Dr_J
        Drr_I
        Drr_J
        It
        Jt
        ArowChange % for storing the row changes made to A and B for imposing BCs
        BrowChange
        rowPlacement % for putting the rows in the right places
    end
    methods
        function obj = scattResComp2d(Nt,Nrs,Rs)
            if mod(Nt,2) == 1
                Nt = Nt + 1; 
                fprintf('Nt must be even. Changing Nt to %d\n', Nt);
            end
            obj.Nt = Nt; obj.Nrs = Nrs; obj.Rs = Rs; obj.Nr = sum(Nrs);
            obj.BCrows = (1:obj.Nt)*obj.Nr;
            
            % mesh between 0 and first breakpoint is half-Cheb type
            [D,x] = cheb(2*Nrs(1)-1,-Rs(1),Rs(1));
            r = x(Nrs(1)+1:end);
            D2 = D*D;
            % blocks we need
            D3  = D2(Nrs(1)+1:end,1:Nrs(1));
            D3t = D2(Nrs(1)+1:end,Nrs(1):-1:1);
            D4  = D2(Nrs(1)+1:end,Nrs(1)+1:end);
            E3  =  D(Nrs(1)+1:end,1:Nrs(1));
            E3t =  D(Nrs(1)+1:end,Nrs(1):-1:1);
            E4  =  D(Nrs(1)+1:end,Nrs(1)+1:end);
            obj.r = r; 
            
            It = eye(Nt); Ithalf = eye(Nt/2);
            Jt = [0*Ithalf, Ithalf; Ithalf, 0*Ithalf];

            obj.Drs_even = E4 + E3t;
            obj.Drs_odd  = E4 - E3t;

            Dr_I  = E4; Dr_J  = E3t; % append onto these
            Drr_I = D4; Drr_J = D3t;
            
            % unscaled version
            [D,x] = cheb(2*Nrs(1)-1);
            E3t =  D(Nrs(1)+1:end,Nrs(1):-1:1);
            E4  =  D(Nrs(1)+1:end,Nrs(1)+1:end);
            obj.Drs_even_unscaled = E4 + E3t;
            obj.Drs_odd_unscaled  = E4 - E3t;
            
            % conversion from cheb coeffs to fourier values
            N = Nrs(1);
            Winv_even = cos( (pi/(2*N-1)) * (-N+1:0).' * (0:2:2*N-2) );
            Winv_odd  = cos( (pi/(2*N-1)) * (-N+1:0).' * (1:2:2*N-1) );
            
            % for the rest of the break points, use usual Cheb mesh
            for ii = 2:length(Rs)
                [Dr,r] = cheb(Nrs(ii)-1,Rs(ii-1),Rs(ii)); 
                obj.Drs_even = blkdiag(obj.Drs_even, Dr);
                obj.Drs_odd  = blkdiag(obj.Drs_odd , Dr);
                obj.r = [obj.r; r];
                Dr_I = blkdiag(Dr_I,  Dr);
                Dr_J = blkdiag(Dr_J,0*Dr);
                
                Drr_I = blkdiag(Drr_I,Dr^2);
                Drr_J = blkdiag(Drr_J,0*Dr);  
                
                % unscaled version
                [Dr,r] = cheb(Nrs(ii)-1); 
                obj.Drs_even_unscaled = blkdiag(obj.Drs_even_unscaled,Dr);
                obj.Drs_odd_unscaled  = blkdiag(obj.Drs_odd_unscaled ,Dr);
                
                N = Nrs(ii) - 1;
                Winv_ii = cos( (pi/N) * (-N:0).' * (0:N) );
                Winv_even = blkdiag(Winv_even, Winv_ii);
                Winv_odd  = blkdiag(Winv_odd , Winv_ii);
            end
            
            % save stuff for fast matvec
            obj.Dr_I = Dr_I;
            obj.Dr_J = Dr_J;
            obj.Drr_I = Drr_I;
            obj.Drr_J = Drr_J;
            obj.It = It;
            obj.Jt = Jt;

            % Maps from values to derivatives on whole mesh
            obj.Dr  = kron(It, Dr_I) + kron(Jt, Dr_J);
            obj.Drr = kron(It,Drr_I) + kron(Jt,Drr_J);
            
            % make Dtt (pulled directly from Spectral Methods in Matlab)
            dtheta = 2*pi/Nt;
            obj.theta = (0:Nt-1)*2*pi/Nt;
            e = ones(Nt,1);
            obj.Dtt = toeplitz([-pi^2/(3*dtheta^2)-1/6 ...
                 .5*(-1).^(2:Nt)./sin(dtheta*(1:Nt-1)/2).^2]);

            [rr,tt] = meshgrid(obj.r, (0:Nt)*2*pi/Nt); % added theta = 2*pi for 3d plotting
            
            [xx,yy] = pol2cart(tt,rr);
            obj.rr = rr; obj.tt = tt; obj.xx = xx; obj.yy = yy;
            
            xs = obj.r*cos(obj.theta); obj.xs = xs(:);
            ys = obj.r*sin(obj.theta); obj.ys = ys(:);
            
            rs = repmat(obj.r,1,obj.Nt);     obj.rs = rs(:);
            ts = repmat(obj.theta,obj.Nr,1); obj.ts = ts(:);
            
            zs = obj.r*exp(1i*obj.theta); obj.zs = zs(:);

            obj.Ns = -obj.Nt/2+1:obj.Nt/2;
            obj.U    = exp( 1i*obj.theta'*obj.Ns);
            obj.Uinv = exp(-1i*obj.Ns'*obj.theta)/obj.Nt;
            obj.Winv_even = Winv_even;
            obj.Winv_odd  = Winv_odd;
            
            % derivatives in Chebyshev basis (simplistic view)
            obj.Drs_cc_even = obj.radialVals2chebCoeffs(obj.Drs_even*obj.Winv_even,'odd');
            obj.Drs_cc_odd  = obj.radialVals2chebCoeffs(obj.Drs_odd *obj.Winv_odd ,'even' );
%             obj.Dr_cc  = obj.chebOperatorFromValuesOperator(obj.Dr);
            obj.Dr_cc = obj.apply_Winv_on_right(obj.Dr*kron(obj.U,eye(obj.Nr)));

            % put interface conditions at edges of subregions
%             obj.valRows = cumsum(obj.Nrs(1:end-1));
%             obj.derRows = obj.valRows + 1;
            obj.derRows = cumsum(obj.Nrs(1:end-1));
            obj.valRows = obj.derRows + 1;
            
            % put interface conditions where coefficients are smallest
            obj.valRows_cc = cumsum(obj.Nrs(1:end-1));
            obj.derRows_cc = obj.valRows_cc + obj.Nrs(2:end) - 1; % right before a valRow (or the BCrow)
            
            % grab negative of old IC/BC rows in A - k^2*B
            globValRows = repmat(obj.valRows',1,obj.Nt) + obj.Nr*ones(length(obj.valRows),1)*(0:obj.Nt-1);
            globDerRows = repmat(obj.derRows',1,obj.Nt) + obj.Nr*ones(length(obj.derRows),1)*(0:obj.Nt-1);
            ICrows = [globValRows; globDerRows]; ICrows = ICrows(:);
            ICBCrows = [reshape(ICrows, 2*length(obj.valRows), obj.Nt); obj.BCrows];
            obj.ICBCrows = ICBCrows(:)';
            
            numFirstIndices = 2*length(obj.valRows) + 1;
            firstBlock = zeros(obj.Nr,numFirstIndices);
            firstIndices = obj.ICBCrows(1:numFirstIndices);
            firstBlock(firstIndices,:) = eye(numFirstIndices);
            obj.rowPlacement = kron(eye(obj.Nt),firstBlock);
        end

        % solve the scattering problem
        scattValuesVec  = solve(obj,k,incfun)
        scattFourierVec = solve_fc(obj,k,incfun)
        scattChebVec    = solve_cc(obj,k,incfun)
        
        % enforce the differential equation in values basis
        enforceDE(obj,Lr,Vs,coords,r,s)        
        enforceIC(obj)

        % enforce the differential equation in fourier basis
        enforceDE_fc(obj)
        enforceIC_fc(obj)
        
        % enforce in Cheb basis
        enforceDE_cc(obj)
        enforceIC_cc(obj)
        
        % stuff for transforming between valuesVec and fourierVec quickly
        y = apply_U(obj,x)
        y = apply_Uinv(obj,x)
        y = apply_U_kron_I(obj,x)
        y = apply_Uinv_kron_I(obj,x)
        y = apply_UT_kron_I(obj,x)
        Y = apply_U_kron_I_on_right(obj,X)
        
        % converting between types of vectors
        valuesVec  = valuesVecFromFourierVec(obj,fourierVec)
        fourierVec = fourierVecFromValuesVec(obj,valuesVec)
        fourierVec = fourierVecFromChebVec(obj,chebVec,deriv)
        chebVec    = chebVecFromFourierVec(obj,fourierVec,deriv)
        valuesVec  = valuesVecFromChebVec(obj,chebVec,deriv)
        chebVec    = chebVecFromValuesVec(obj,valuesVec,deriv)
        
        % helpers for cheb basis
        c = radialVals2chebCoeffs(obj,v,parity)
        v = chebCoeffs2radialVals(obj,c,parity)
        X = apply_W(obj,X,deriv)
        X = apply_Winv_on_right(obj,X)
        X = chebOperatorFromValuesOperator(obj,X)

        X = fourierOperatorFromValuesOperator(obj,X)
        
        % vectors from functions
        valuesVec = valuesVecFromFunCellArray(obj,fs,coords)  
        fourierVec = fourierVecFromFunCellArray(obj,fs,coords)
        valuesVec = valuesVecFromFun(obj,f,coords)
        fourierVec = fourierVecFromFun(obj,f,coords)
        
        % right-hand side for scattering problem
        RHS = RHSfromFun(obj,incfun)   
        RHS = RHSfromFun_cc(obj,incfun)
        RHS = RHSfromFun_fc(obj,incfun)
        
        % plotting vectors
        u = plotFourierVec(obj,fourierVec,titlestr,part)
        u = plotValuesVec(obj,valuesVec,titlestr,part)
        
        % plot function values along a ray
        plotRayValuesVec(obj,rayValuesVec,style,name)
        
        % plotting functions
        u = plotFun(obj,f,titlestr,part,coords)
        
        % eigenvalue count stuff
        dfoverf = dLogTz(obj,z,br)
    end
    methods (Static)
        % helpers for cheb basis
        c = halfChebVals2chebCoeffs(v,parity)
        c = fullChebVals2chebCoeffs(v)
        v = chebCoeffs2halfChebVals(c,parity)
        v = chebCoeffs2fullChebVals(c)
        
        function [evecs,evals] = sort_eigs(evecs,evals)
            % sorts from left to right
            evals = diag(evals);
            [~,indx] = sort(real(evals));
            evals = evals(indx);
            evecs = evecs(:,indx);            
        end
        
        function n_closest_p = get_closest(p,n,p0)
            % gets n points in p closest to p0
            distances = abs(p-p0);
            [~,idx] = sort(distances);
            sorted_p = p(idx);
            n_closest_p = sorted_p(1:n);
        end
        
        % bessel functions and derivatives with respect to x
        function Jnkx = J(n,k,x)
            Jnkx = besselj(n,k*x);   
        end
        function dJnkx = dJ(n,k,x)
            dJnkx = k*(scattResComp2d.J(n-1,k,x)-scattResComp2d.J(n+1,k,x))/2;
        end
        function Ynkx = Y(n,k,x)
            Ynkx = bessely(n,k*x);   
        end
        function dYnkx = dY(n,k,x)
            dYnkx = k*(scattResComp2d.Y(n-1,k,x)-scattResComp2d.Y(n+1,k,x))/2;
        end
        function Hnkx = H(n,k,x)
            Hnkx = besselh(n,1,k*x);
        end
        function dHnkx = dH(n,k,x)
            dHnkx = k*(scattResComp2d.H(n-1,k,x)-scattResComp2d.H(n+1,k,x))/2;
        end
        function y = apply_kron(X,Y,xhat)
            % y = kron(X,Y)*xhat(:)
            y = Y*xhat*X.';  y = y(:);
        end
    end
end

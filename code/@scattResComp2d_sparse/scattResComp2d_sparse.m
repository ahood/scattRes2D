classdef scattResComp2d_sparse < handle
    % General framework for scattering computation and plotting.
    properties
        % Only commenting properties that are not in scattResComp2d. 
        Nt     
        Nrs    
        Nr     
        Rs     
        r      
        rs
        ts
        xs
        ys
        zs
        rr
        tt
        xx
        yy
        Dtt    
        Dtt_star
        theta  
        BCrows
        valRows
        derRows
        globalICBCrows % the locations of all rows changed for ICs or BCs
        localICBCrows  % rows changed for ICs or BCs in a given block
        Ns
        E4  % kron(I,E4) + kron(J,E3t) maps values to r derivatives on half cheb mesh
        D4  % kron(I,D4) + kron(J,D3t) maps values to r second derivatives on half cheb mesh
        E3t % above
        D3t % above
        E4_star
        D4_star
        E3t_star
        D3t_star
        Dr_I
        Drr_I
        Dr_J
        Drr_J
        R
        U
        Uinv
        Winv_even
        Winv_odd
        ArowChange % for storing the row changes made to A and B for imposing BCs
        BrowChange
        localPlacement % for putting the rows in the right places
    end
    methods
        function obj = scattResComp2d_sparse(Nt,Nrs,Rs)
            if mod(Nt,2) == 1
                Nt = Nt + 1; 
                fprintf('Nt must be even. Changing Nt to %d\n', Nt);
            end
            obj.Nt = Nt; obj.Nrs = Nrs; obj.Rs = Rs; obj.Nr = sum(Nrs);
            obj.BCrows = (1:obj.Nt)*obj.Nr;
%             obj.valRows = cumsum(obj.Nrs(1:end-1));
%             obj.derRows = obj.valRows + 1;
            obj.derRows = cumsum(obj.Nrs(1:end-1));
            obj.valRows = obj.derRows + 1;
            
            % mesh between 0 and first breakpoint is half-Cheb type
            [D,x] = cheb(2*Nrs(1)-1,-Rs(1),Rs(1));
            r = x(Nrs(1)+1:end);
            D2 = D*D;
            % blocks we need
            obj.D3t = D2(Nrs(1)+1:end,Nrs(1):-1:1);
            obj.D4  = D2(Nrs(1)+1:end,Nrs(1)+1:end);
            obj.E3t =  D(Nrs(1)+1:end,Nrs(1):-1:1);
            obj.E4  =  D(Nrs(1)+1:end,Nrs(1)+1:end);
            obj.D3t_star = obj.D3t';
            obj.D4_star  = obj.D4';
            obj.E3t_star = obj.E3t';
            obj.E4_star  = obj.E4';

            obj.r = r; 
            Dr_I = obj.E4;  Drr_I = obj.D4; % append onto these
            Dr_J = obj.E3t; Drr_J = obj.D3t;
            
            N = Nrs(1);
            Winv_even = cos( (pi/(2*N-1)) * (-N+1:0).' * (0:2:2*N-2) );
            Winv_odd  = cos( (pi/(2*N-1)) * (-N+1:0).' * (1:2:2*N-1) );
            
            % for the rest of the break points, use usual Cheb mesh
            for ii = 2:length(Rs)
                [Dr,r] = cheb(Nrs(ii)-1,Rs(ii-1),Rs(ii)); 
                obj.r = [obj.r; r];
                Dr_I  = blkdiag(Dr_I ,Dr  );
                Drr_I = blkdiag(Drr_I,Dr^2);
                Dr_J = blkdiag(Dr_J,0*Dr);
                Drr_J = blkdiag(Drr_J,0*Dr); 
                
                N = Nrs(ii) - 1;
                Winv_ii = cos( (pi/N) * (-N:0).' * (0:N) );
                Winv_even = blkdiag(Winv_even, Winv_ii);
                Winv_odd  = blkdiag(Winv_odd , Winv_ii);
            end
            obj.Dr_I = Dr_I; obj.Drr_I = Drr_I;    
            obj.Dr_J = Dr_J; obj.Drr_J = Drr_J;

            obj.R = spdiags(1./obj.r,0,obj.Nr,obj.Nr);

            dtheta = 2*pi/Nt;  
            obj.theta = (0:Nt-1)*2*pi/Nt;
            e = ones(Nt,1);
            obj.Dtt = toeplitz([-pi^2/(3*dtheta^2)-1/6 ...
                 .5*(-1).^(2:Nt)./sin(dtheta*(1:Nt-1)/2).^2]);
            obj.Dtt_star = obj.Dtt';
                        
            obj.Ns = -obj.Nt/2+1:obj.Nt/2;
            
            [rr,tt] = meshgrid(obj.r,(0:Nt)*2*pi/Nt);
            
            [xx,yy] = pol2cart(tt,rr);
            obj.rr = rr; obj.tt = tt; obj.xx = xx; obj.yy = yy;            
            
            xs = obj.r*cos(obj.theta); obj.xs = xs(:);
            ys = obj.r*sin(obj.theta); obj.ys = ys(:);            

            rs = repmat(obj.r,1,obj.Nt);     obj.rs = rs(:); 
            ts = repmat(obj.theta,obj.Nr,1); obj.ts = ts(:);
            obj.zs = obj.xs + 1i*obj.ys;
            
            obj.U    = exp( 1i*obj.theta'*obj.Ns);
            obj.Uinv = exp(-1i*obj.Ns'*obj.theta)/obj.Nt;
            obj.Winv_even = Winv_even;
            obj.Winv_odd  = Winv_odd;
            
            % locations of ICs and BCs
            globValRows = repmat(obj.valRows',1,obj.Nt) + obj.Nr*ones(length(obj.valRows),1)*(0:obj.Nt-1);
            globDerRows = repmat(obj.derRows',1,obj.Nt) + obj.Nr*ones(length(obj.derRows),1)*(0:obj.Nt-1);
            ICrows = [globValRows; globDerRows]; ICrows = ICrows(:);
            ICBCrows = [reshape(ICrows, 2*length(obj.valRows), obj.Nt); obj.BCrows];
%             obj.globalICBCrows = ICBCrows(:)';
            obj.globalICBCrows = sort(ICBCrows(:))';
            
%             localRows = [obj.valRows; obj.derRows]; localRows = localRows(:);
            localRows = [obj.valRows; obj.derRows]; localRows = sort(localRows(:));
            localRows(end+1) = obj.Nr;
%             obj.localICBCrows = localRows(:)';
            obj.localICBCrows = sort(localRows(:))';
            
            localPlacement = spalloc(obj.Nr,2*length(obj.valRows) + 1, length(localRows));
            localPlacement(localRows,:) = speye(length(localRows));
            obj.localPlacement = localPlacement;
        end
        
        % solve the scattering problem
        scattValuesVec = solve(obj,k,incfun,precond,verbose)
        
        function y = apply_DtNmap(obj,x,k)
            % y = U*D*U\x for D diagonal matrix of DtN map coefficients
            y = obj.apply_Uinv(x);
            y = obj.DtNcoeffs(obj.Ns',k,obj.r(end)).*y;
            y = obj.apply_U(y);
        end
        
        function yhat = apply_Drr_J(obj,xhat)
            % yhat = Drr_J*xhat where xhat could be a matrix
            yhat = [obj.D3t*xhat(1:obj.Nrs(1)    ,:); ...
                          0*xhat(obj.Nrs(1)+1:end,:)];
        end
    
        function yhat = apply_Drr_J_star(obj,xhat)
            % yhat = Drr_J'*xhat where xhat could be a matrix
            yhat = [obj.D3t_star*xhat(1:obj.Nrs(1)    ,:); ...
                               0*xhat(obj.Nrs(1)+1:end,:)];
        end
    
        function y = apply_J_kron_Drr_J(obj,xhat)
            % y = kron(J,Drr_J)*x, x = xhat(:)
            
            xhat = obj.apply_Drr_J(xhat);

            cols1 = 1:obj.Nt/2;
            xhat = [xhat(:,cols1(end)+1:end), xhat(:,cols1)];
            
            y = xhat(:);
        end 
        
        function y = apply_J_kron_Drr_J_star(obj,xhat)
            % y = kron(J,Drr_J')*x, x = xhat(:)
            
            xhat = obj.apply_Drr_J_star(xhat);

            cols1 = 1:obj.Nt/2;
            xhat = [xhat(:,cols1(end)+1:end), xhat(:,cols1)];
            
            y = xhat(:);
        end 
        
        function yhat = apply_Dr_J(obj,xhat)
            % yhat = Dr_J*xhat where xhat could be a matrix
            
            yhat = [obj.E3t*xhat(1:obj.Nrs(1)    ,:); ...
                          0*xhat(obj.Nrs(1)+1:end,:)];
        end
        
        function yhat = apply_Dr_J_star(obj,xhat)
            % yhat = Dr_J'*xhat where xhat could be a matrix
            
            yhat = [obj.E3t_star*xhat(1:obj.Nrs(1)    ,:); ...
                               0*xhat(obj.Nrs(1)+1:end,:)];
        end
        
        function y = apply_J_kron_Dr_J(obj,xhat)
            % y = kron(J,Dr_J)*x, x = xhat(:)
            
            xhat = obj.apply_Dr_J(xhat);

            cols1 = 1:obj.Nt/2;
            xhat = [xhat(:,cols1(end)+1:end), xhat(:,cols1)];
            
            y = xhat(:);
        end
        
        function y = apply_J_kron_Dr_J_star(obj,xhat)
            % y = kron(J,Dr_J')*x, x = xhat(:)
            
            xhat = obj.apply_Dr_J_star(xhat);

            cols1 = 1:obj.Nt/2;
            xhat = [xhat(:,cols1(end)+1:end), xhat(:,cols1)];
            
            y = xhat(:);
        end
        
        function xhat = apply_Dr_I(obj,xhat)
            % xhat = Dr_I*xhat where xhat could be a matrix
            
            % apply the half cheb mesh part
            rows = 1:obj.Nrs(1);
            xhat(rows,:) = obj.E4*xhat(rows,:);
            
            % apply the regular cheb mesh parts
            for n = 2:length(obj.Nrs)
                rows = rows(end) + (1:obj.Nrs(n));
                a = obj.Rs(n-1); b = obj.Rs(n);
                xhat(rows,:) = chebfft_ah(xhat(rows,:),a,b);
            end
        end
        
        function xhat = apply_Dr_I_star(obj,xhat)
            % xhat = Dr_I'*xhat where xhat could be a matrix
            
            % apply the half cheb mesh part
            rows = 1:obj.Nrs(1);
            xhat(rows,:) = obj.E4_star*xhat(rows,:);
            
            % apply the regular cheb mesh parts
            for n = 2:length(obj.Nrs)
                rows = rows(end) + (1:obj.Nrs(n));
                a = obj.Rs(n-1); b = obj.Rs(n);
                xhat(rows,:) = chebfft_star_ah(xhat(rows,:),a,b);
            end
        end
        
        function y = apply_I_kron_Dr_I(obj,xhat)
            % y = kron(I,Dr_I)*x, x = xhat(:)
            
            xhat = obj.apply_Dr_I(xhat);
            y = xhat(:);
        end
        
        function y = apply_I_kron_Dr_I_star(obj,xhat)
            % y = kron(I,Dr_I')*x, x = xhat(:)
            
            xhat = obj.apply_Dr_I_star(xhat);
            y = xhat(:);
        end        
        
        function xhat = apply_Drr_I(obj,xhat)
            % xhat = Drr_I*xhat where xhat could be a matrix
 
            % apply the half cheb mesh part
            rows = 1:obj.Nrs(1);
            xhat(rows,:) = obj.D4*xhat(rows,:);
            
            for n = 2:length(obj.Nrs)
                rows = rows(end) + (1:obj.Nrs(n)); 
                a = obj.Rs(n-1); b = obj.Rs(n);
                xhat(rows,:) = chebfft_ah(chebfft_ah(xhat(rows,:),a,b),a,b);
            end            
        end

        function xhat = apply_Drr_I_star(obj,xhat)
            % xhat = Drr_I'*xhat where xhat could be a matrix
 
            % apply the half cheb mesh part
            rows = 1:obj.Nrs(1);
            xhat(rows,:) = obj.D4_star*xhat(rows,:);
            
            for n = 2:length(obj.Nrs)
                rows = rows(end) + (1:obj.Nrs(n)); 
                a = obj.Rs(n-1); b = obj.Rs(n);
                xhat(rows,:) = chebfft_star_ah(chebfft_star_ah(xhat(rows,:),a,b),a,b);
            end            
        end

        function y = apply_I_kron_Drr_I(obj,xhat)
            % y = kron(I,Drr_I)*x, x = xhat(:)
            
            xhat = obj.apply_Drr_I(xhat);
            y = xhat(:);
        end        
        
        function y = apply_I_kron_Drr_I_star(obj,xhat)
            % y = kron(I,Drr_I')*x, x = xhat(:)
            
            xhat = obj.apply_Drr_I_star(xhat);
            y = xhat(:);
        end        
        
        % do not form, but apply
        [x,DrIxhat,DrJxhat] = apply_op_no_bcs(obj,x0,E0,VvaluesVec)
        x = apply_Aorig_star(obj,x0,conjVvaluesVec)
        
        % have to form boundary condition adjustments -- call once
        % VvaluesVec set up
        ICBCmatrices(obj)
                
        % stuff for transforming between valuesVec and fourierVec quickly
        y = apply_U(obj,x)
        y = apply_Uinv(obj,x)
        y = apply_U_kron_I(obj,x)
        y = apply_Uinv_kron_I(obj,x)
        
        % converting between types of vectors
        valuesVec = valuesVecFromFourierVec(obj,fourierVec)
        fourierVec = fourierVecFromValuesVec(obj,valuesVec)
        fourierVec = fourierVecFromChebVec(obj,chebVec)
        chebVec = chebVecFromFourierVec(obj,fourierVec)
        valuesVec = valuesVecFromChebVec(obj,chebVec)
        chebVec = chebVecFromValuesVec(obj,valuesVec)
        
        % helpers for cheb basis
        c = radialVals2chebCoeffs(obj,v,parity)
        v = chebCoeffs2radialVals(obj,c,parity)
        
        % vectors from functions
        valuesVec = valuesVecFromFunCellArray(obj,fs,coords)  
        fourierVec = fourierVecFromFunCellArray(obj,fs,coords)
        valuesVec = valuesVecFromFun(obj,f,coords)
        fourierVec = fourierVecFromFun(obj,f,coords)
        
        % right-hand side for scattering problem
        RHS = RHSfromFun(obj,incfun)   
        
        % plotting vectors
        u = plotFourierVec(obj,fourierVec,titlestr,part)
        u = plotValuesVec(obj,valuesVec,titlestr,part)
        
        % plot function values along a ray
        f = plotRayValuesVec(obj,rayValuesVec,style,name)
        
        % plotting functions
        u = plotFun(obj,f,titlestr,part,coords)
    end
    methods (Static)
        function [evecs,evals] = sort_eigs(evecs,evals)
            evals = diag(evals);
            [~,indx] = sort(real(evals));
            evals = evals(indx);
            evecs = evecs(:,indx);            
        end
        
        % helpers for cheb basis
        c = halfChebVals2chebCoeffs(v,parity)
        c = fullChebVals2chebCoeffs(v)
        v = chebCoeffs2halfChebVals(c,parity)
        v = chebCoeffs2fullChebVals(c)
        
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

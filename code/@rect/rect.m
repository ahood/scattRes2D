classdef rect < handle
    properties
        x1
        x2
        y1
        y2
        nx
        ny
        X
        Y
        Z
    end
    methods
        function obj = rect(x1,x2,y1,y2,nx,ny)
            obj.x1 = x1; obj.x2 = x2; obj.y1 = y1; obj.y2 = y2;
            if (~isempty(nx) && ~isempty(ny))
                obj.set_mesh(nx,ny); 
                obj.nx = nx; obj.ny = ny;
            end
        end
        function set_mesh(obj,nx,ny)
            [X,Y] = meshgrid( linspace(obj.x1,obj.x2,nx), ...
                              linspace(obj.y1,obj.y2,ny) );
            obj.X = X; obj.Y = Y; obj.Z = X + 1i*Y;
        end
        function tf = contains(obj,pts)
            res = real(pts); ims = imag(pts);
            tf = (obj.x1 <= res) & (res <= obj.x2) & ...
                 (obj.y1 <= ims) & (ims <= obj.y2);
        end
        function focus(obj,buff)
            if nargin == 1, buff = 0; end
            axis([obj.x1-buff obj.x2+buff obj.y1-buff obj.y2+buff])
        end
        function log10contour(obj,f,newfig,labels,v,color,verbose)
            if nargin < 3 || isempty(newfig), newfig = 1; end
            if nargin < 4 || isempty(labels), labels = 0; end
            if nargin < 6 || isempty(color), color = 'k'; end
            if nargin < 7 || isempty(verbose), verbose = 0; end
            fvals = 0*obj.Z;
            if verbose
                [m,n] = size(obj.Z); N = m*n;
                for ii = 1:N
                    if mod(ii,m) == 0, 
                        fprintf('Did %d out of %d points (%.2f done)\n', ii, N, ii/N);
                    end
                    z = obj.Z(ii);
                    fvals(ii) = f(z);
                end
            else
                for ii = 1:numel(obj.Z)
                    z = obj.Z(ii);
                    fvals(ii) = f(z);
                end
            end    
            
            if newfig, figure; end
            if nargin < 5 || isempty(v)
                [C,h] = contour(obj.X,obj.Y,log10(fvals));
            else
                if v(1) == v(2)
                    [C,h] = contour(obj.X,obj.Y,log10(fvals),v,'linecolor',color,'linewidth',2);
                    return
                else
                    [C,h] = contour(obj.X,obj.Y,log10(fvals),v);
                end
            end
            if labels, clabel(C,h,'fontsize',12); else colorbar; end
        end
        function contour(obj,f)
            fvals = 0*obj.Z;
            parfor ii = 1:numel(obj.Z)
                z = obj.Z(ii);
                fvals(ii) = f(z);
            end
            figure, contour(obj.X,obj.Y,fvals); colorbar;
        end
        function draw(obj)
            plot([obj.x1 obj.x1 obj.x2 obj.x2 obj.x1], ...
                 [obj.y1 obj.y2 obj.y2 obj.y1 obj.y1]);
        end
    end
end
classdef ratApprox < handle
    properties
        region % region on which we make rat approx
        f
        z % nodes
        w % weights
    end
    methods
        function obj = ratApprox(region,f,t)
            % get rat approx to f where t is mesh of [0,1]
            if nargin == 3
                if isempty(t), t = linspace(0,1,100); end
                obj.region = region;
                obj.f = f;
                obj.z = region.phi(t);
                h = [t(2)-t(1), t(3:end)-t(1:end-2), t(end)-t(end-1)];
                obj.w = f(obj.z).*(region.dphi(t)).*h/4/pi/1i;
            end
        end
        function fzapprox = eval_at(obj,z0)
            fzapprox = sum( obj.w./(obj.z - z0) );
        end
        function E = compute_error(obj,r)
            % computes error on rectangle with obj rep r
            E = 0*r.Z;
            for ii = 1:numel(r.Z)
                z = r.Z(ii);
                E(ii) = abs( obj.f(z) - obj.eval_at(z) );
            end
        end
        function show_error(obj,r,newfig,labels,v,color,verbose)
            if nargin < 3, newfig = 1; end
            if nargin < 4, labels = 0; end
            if nargin < 5, v = []; end
            if nargin < 6, color = 'k'; end
            if nargin < 7, verbose = 0; end
            r.log10contour(@(z) abs(obj.f(z) - obj.eval_at(z) ), newfig, labels, v, color,verbose);
            title('absolute error');
        end
    end
end
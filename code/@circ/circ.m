classdef circ < handle
    properties
        c % center
        r % radius
        X
        Y
        Z
        bb % bounding box
    end
    methods
        function obj = circ(c,r,nx,ny)
            obj.c = c; obj.r = r;
            obj.bb = rect(real(c)-r,real(c)+r,...
                          imag(c)-r,imag(c)+r,...
                          [],[]);
            if (~isempty(nx) && ~isempty(ny)), obj.set_mesh(nx,ny); end
        end
        function set_mesh(obj,nx,ny)
            c = obj.c; r = obj.r;
            x1 = real(c)-r; x2 = real(c)+r;
            y1 = imag(c)-r; y2 = imag(c)+r;
            [X,Y] = meshgrid( linspace(x1,x2,nx), linspace(y1,y2,ny) );
            obj.X = X; obj.Y = Y; obj.Z = X + 1i*Y;
        end
        function tf = contains(obj,pts)
            tf = abs( pts - obj.c ) < r;
        end
        function phit = phi(obj,t)
            phit = obj.c + obj.r*exp(2*pi*1i*t);
        end
        function dphit = dphi(obj,t)
            dphit = obj.r*2*pi*1i*exp(2*pi*1i*t);
        end
        function draw(obj)
            t = linspace(0,1,2*pi*obj.r*1e2);
            p = obj.phi(t);
            plot(p)
        end
        function focus(obj,buff)
            obj.bb.focus(buff)
        end
    end
end
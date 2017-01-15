classdef ellipse < handle
    properties
        c % center
        theta % rotation angle
        a % axis length along theta direction
        b % other one
        f1 % one focus
        f2 % the other one
        R % |z - f1| + |z - f2| = R
        X
        Y
        Z
        bb % bounding box
    end
    methods
        function obj = ellipse(c,theta,a,b,nx,ny)
            obj.c = c; obj.theta = theta; obj.a = a; obj.b = b;
            obj.bb = rect(real(c)-a,real(c)+a,...
                          imag(c)-b,imag(c)+b,...
                          nx,ny);
            f = sqrt(max(a,b)^2 - min(a,b)^2);
            if a > b
                obj.f1 =  f*exp(1i*theta) + c;
                obj.f2 = -f*exp(1i*theta) + c;
            else
                obj.f1 =  1i*f*exp(1i*theta) + c;
                obj.f2 = -1i*f*exp(1i*theta) + c;
            end
            obj.R = 2*max(a,b);
            if (~isempty(nx) && ~isempty(ny)), obj.set_mesh(nx,ny); end
        end
        function set_mesh(obj,nx,ny)
            c = obj.c; r = max(obj.a,obj.b);
            x1 = real(c)-r; x2 = real(c)+r;
            y1 = imag(c)-r; y2 = imag(c)+r;
            [X,Y] = meshgrid( linspace(x1,x2,nx), linspace(y1,y2,ny) );
            obj.X = X; obj.Y = Y; obj.Z = X + 1i*Y;
        end
        function tf = contains(obj,pts)
            tf = abs(pts - obj.f1) + abs(pts - obj.f2) < obj.R;
        end
        function phit = phi(obj,t)
            phit = obj.c + exp(1i*obj.theta)*(obj.a*cos(2*pi*t) + 1i*obj.b*sin(2*pi*t));
        end
        function dphit = dphi(obj,t)
            dphit = 2*pi*exp(1i*obj.theta)*(-obj.a*sin(2*pi*t) + 1i*obj.b*cos(2*pi*t));
        end
        function draw(obj)
            % ds = r*dtheta with ds < 1e-2 means
            % dtheta = 1e-2/a. Since dtheta = 2*pi/N,
            % N = 2pi/dtheta, or total number of
            % thetas is N = 2*pi*a*1e2
            t = linspace(0,1,2*pi*obj.a*1e2);
            p = obj.phi(t);
            plot(p)
        end
        function focus(obj,buff)
            obj.bb.focus(buff)
        end
    end
end
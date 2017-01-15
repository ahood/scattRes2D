classdef ratApproxWithPoles < ratApprox
    % difference here is that we subtract away contour integrals around
    % poles
    properties
        poles;
        rs; % radii around poles
    end
    methods
        function obj = ratApproxWithPoles(region,f,poles,rs,t)
            obj@ratApprox(region,f,t);
            obj.poles = poles;
            obj.rs = rs;
            % mesh of boundary centered at pole - at least size 20, and
            %   otherwise 
            Ns = max(40,round(100*rs) );
            for jj = 1:length(poles)
                t = linspace(0,1,Ns(jj));
                pDisk = circ(poles(jj), rs(jj), [],[]);
                pRatApprox = ratApprox(pDisk,f,t);
                obj.z = [obj.z,  pRatApprox.z];
                obj.w = [obj.w, -pRatApprox.w];
            end
        end     
    end
end
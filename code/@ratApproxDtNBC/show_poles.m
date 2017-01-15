function show_poles(obj,r)
    colors = jet(obj.dtn.Nt);
%     for n = 0:obj.dtn.Ns
    for n = 0:obj.dtn.Nt/2
        polesn = ratApproxDtNBC.get_poles_in_range(n,r,obj.dtn.r(end),obj.br);
        if obj.br == 1  % uses first branch of sqrt
            resid = norm( besselh(n,1, sqrt(polesn)*obj.dtn.r(end)) );  
        else
            resid = norm( besselh(n,1,-sqrt(polesn)*obj.dtn.r(end)) );
        end
%                 fprintf('for n = %d, residual wrt H is %4.2e\n', n, resid);
        plot(polesn,'.','color',colors(n+1,:),'markersize',14,'displayname',num2str(n));
    end
    axis([r.x1 r.x2 r.y1 r.y2])
    legend show
end
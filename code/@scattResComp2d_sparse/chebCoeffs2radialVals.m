function v = chebCoeffs2radialVals(obj,c,parity)
% See radialVals2chebCoeffs. This does opposite.

    v = [];
    
    j = 1;
    cj = c(1:obj.Nrs(j),:); c = c(obj.Nrs(j)+1:end,:);
    vj = scattResComp2d.chebCoeffs2halfChebVals(cj,parity); 
    v = [v; vj];
    
    for j = 2:length(obj.Nrs)
        cj = c(1:obj.Nrs(j),:); c = c(obj.Nrs(j)+1:end,:);
        vj = scattResComp2d.chebCoeffs2fullChebVals(cj); 
        v = [v; vj];
    end

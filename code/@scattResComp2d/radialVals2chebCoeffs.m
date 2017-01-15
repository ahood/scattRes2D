function c = radialVals2chebCoeffs(obj,v,parity)
% Converts v = f_j(obj.r) to c = c_j as described in chebVecFromFourierVec,
% with f_j an even or odd function specified by parity = 'even' or 'odd'.
% The variable v can either represent a vector or a matrix whose columns
% are the f_j(obj.r).

    c = [];
    
    j = 1;
    vj = v(1:obj.Nrs(j),:); v = v(obj.Nrs(j)+1:end,:);
    cj = scattResComp2d.halfChebVals2chebCoeffs(vj,parity); 
    c = [c; cj];
    
    for j = 2:length(obj.Nrs)
        vj = v(1:obj.Nrs(j),:); v = v(obj.Nrs(j)+1:end,:);
        cj = scattResComp2d.fullChebVals2chebCoeffs(vj); 
        c = [c; cj];
    end
    
function v = chebCoeffs2halfChebVals(c,parity)
% See halfChebVals2chebCoeffs. This does opposite.

    N = size(c,1);
    C = zeros(2*N,size(c,2));
    if strcmp(parity,'even')
        C(1:2:end,:) = c;
    elseif strcmp(parity,'odd')
        C(2:2:end,:) = c;
    end
    
    V = scattResComp2d.chebCoeffs2fullChebVals(C);
    
    v = V(N+1:end,:);
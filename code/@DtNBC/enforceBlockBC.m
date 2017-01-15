function X = enforceBlockBC(obj,blk,X)
% takes a diagonal block associated to azimuth n and enforces interface and
% boundary conditions
%     if mod(n,2), Dr = obj.Dr_as; else Dr = obj.Dr_s; end
    
    shift = (blk-1)*obj.Nr;

    % boundary condition
    X(end,:) = obj.Dr(shift + obj.Nr,:);
end

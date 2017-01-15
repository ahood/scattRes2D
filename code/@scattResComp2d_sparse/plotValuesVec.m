function u = plotValuesVec(obj,valuesVec,titlestr,part)
    if nargin < 4, part = @real; end
    
    u = reshape(part(valuesVec),obj.Nr,[]); % one column per ray
    u = u.'; % transpose to use mesh()
    u = [u; u(1,:)]; % so the plotted surface isn't missing a sector

    mesh(obj.xx,obj.yy,u); view(30,40)
    title(titlestr)

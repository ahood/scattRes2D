function f = plotRayValuesVec(obj,rayValuesVec,style,name,part)
% Takes the values of a function on one ray, that is, on the r mesh and at
% a particular theta value, turn it into a chebfun, and plot using chebfun.
    if nargin < 5, part = @real; end
    
    % cell array to hold piecewise-defined parts
    c = cell(1,length(obj.Rs));

    % deal with half cheb mesh part
    y = rayValuesVec(1:obj.Nrs(1));
    y_aug = [flipud(y); y];
    c{1} = y_aug;

    % deal with other parts
    for n = 2:length(obj.Rs)
        y = rayValuesVec(sum(obj.Nrs(1:n-1))+1:end);
        c{n} = y;
    end

    % make a chebfun and plot it
    f = chebfun(c, [-obj.Rs(1), obj.Rs]);
    plot(part(f),style,'displayname',name);
end
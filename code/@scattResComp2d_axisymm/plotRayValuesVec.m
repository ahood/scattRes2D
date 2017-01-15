function plotRayValuesVec(obj,rayValuesVec,style,name)
    if nargin < 4, part = @real; end

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
    plot(f,style,'displayname',name);

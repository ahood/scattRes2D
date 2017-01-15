function valuesVec = valuesVecFromFun(obj,f,coords)
% For f a function of two variables, creates a vector of the form
%    [f(obj.r,obj.theta(1)); f(obj.r,obj.theta(2)); ...; f(obj.r,obj.theta(end))]
% where obj.r and obj.theta are meshes in the r and theta directions,
% respectively.
% The choices for coords are:
% - 'rect', meaning user passes f = f(x,y)
% - 'polar', for f = f(r,theta)
% - 'complex', for f = f(z).

    if     strcmp(coords,'rect'   ), valuesVec = f(obj.xs,obj.ys);
    elseif strcmp(coords,'polar'  ), valuesVec = f(obj.rs,obj.ts);
    elseif strcmp(coords,'complex'), valuesVec = f(obj.zs);
    else error('Coordinate type entered isn''t recognized'); end

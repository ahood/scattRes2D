function valuesVec = valuesVecFromFunCellArray(obj,fs,coords)
% For some function f(r,theta) piecewise defined in the r direction, fs is
% a cell array containing the function definitions on
% annuli [0,r1], [r1,r2], etc. 
% For f given above, valuesVec is a vector of the form
%    [f(obj.r,obj.theta(1)); f(obj.r,obj.theta(2)); ...; f(obj.r,obj.theta(end))]
% where obj.r and obj.theta are meshes in the r and theta directions,
% respectively.
% The choices for coords are:
% - 'rect', meaning user passes f = f(x,y)
% - 'polar', for f = f(r,theta)
% - 'complex', for f = f(z).

    valuesRect = zeros(obj.Nr,obj.Nt);
    cumNrs = cumsum(obj.Nrs);
    for ii = 1:length(fs)
        valuesVecii = obj.valuesVecFromFun(fs{ii},coords);
        valuesRectii = reshape(valuesVecii,obj.Nr,[]); 
        rowsii = cumNrs(ii)-obj.Nrs(ii)+1:cumNrs(ii);
        valuesRect(rowsii,:) = valuesRectii(rowsii,:);
    end
    valuesVec = valuesRect(:);

function fourierVec = fourierVecFromFunCellArray(obj,fs,coords)
% For some function f(r,theta) piecewise defined in the r direction, fs is
% a cell array containing the function definitions on
% annuli [0,r1], [r1,r2], etc. 
% For f given above, fourierVec is a vector of the form
%    [f_{-n+1}(obj.r); f_{-n+2}(obj.r); ... ; f_{n}(obj.r)]
% where n = obj.Nt/2 (half the number of points in the theta mesh), 
% obj.r and obj.theta are meshes in the r and theta
% directions, respectively, and the f_{k} are defined by 
%    f(r,theta) = sum_k f_{k}(rs) exp(1i*k*theta)   (Fourier series).
% The choices for coords are:
% - 'rect', meaning user passes f = f(x,y)
% - 'polar', for f = f(r,theta)
% - 'complex', for f = f(z).

    valuesVec = obj.valuesVecFromFunCellArray(fs,coords);
    fourierVec = obj.fourierVecFromValuesVec(valuesVec);

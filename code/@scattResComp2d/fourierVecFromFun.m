function fourierVec = fourierVecFromFun(obj,f,coords)
% For f a function of two variables, creates a vector of the form
%    [f_{-n}(obj.r); f_{-n+1}(obj.r); ...; f_{n}(obj.r)]
% where obj.r is a mesh in the r direction, n = obj.Nt/2 (half the number 
% of mesh points in the theta direction), and the f_{k} are defined by
%    f(r,theta) = sum_k f_{k}(rs) exp(1i*k*theta) (Fourier series).
% The choices for coords are:
% - 'rect', meaning user passes f = f(x,y)
% - 'polar', for f = f(r,theta)
% - 'complex', for f = f(z).

    valuesVec = obj.valuesVecFromFun(f,coords);
    fourierVec = obj.fourierVecFromValuesVec(valuesVec);

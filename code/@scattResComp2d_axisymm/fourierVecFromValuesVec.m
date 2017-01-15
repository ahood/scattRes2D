function fourierVec = fourierVecFromValuesVec(obj,valuesVec)
% Takes a vector of the form 
%    [f(obj.r,obj.theta(1)); f(obj.r,obj.theta(2)); ...; f(obj.r,obj.theta(end))] 
% and transforms to 
%    [f_{-n+1}(obj.r); f_{-n+2}(obj.r); ... ; f_{n}(obj.r)]
% where n = N/2, obj.r and obj.theta are meshes in the r and theta
% direction, respectively, and the f_{k} are defined by 
%    f(r,theta) = sum_k f_{k}(rs) exp(1i*k*theta)   (Fourier series).

    valuesRect = reshape(valuesVec,obj.Nr,obj.Nt);
    fourierVec = obj.apply_Uinv_kron_I(valuesRect);

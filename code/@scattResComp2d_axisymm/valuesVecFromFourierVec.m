function valuesVec = valuesVecFromFourierVec(obj,fourierVec)
% Takes a vector of the form 
%    [f_{-maxn+1}(obj.r); f_{-maxn+2}(obj.r); ... ; f_{maxn}(obj.r)]
% and transforms to 
%    [f(obj.r,obj.theta(1)); f(obj.r,obj.theta(2)); ...; f(obj.r,obj.theta(end))]   
% where obj.r and obj.theta are meshes in the r and theta directions,
% respectively, and the f_{k} are defined by
%    f(r,theta) = sum_k f_{k}(rs) exp(1i*k*theta)  (Fourier series).

    % first figure out how many coefficients there are
    nCoeffs = length(fourierVec)/obj.Nr;
    fourierRect = reshape(fourierVec,obj.Nr,[]);    
    if nCoeffs == obj.Nt
        valuesVec = obj.apply_U_kron_I(fourierRect);
    else
        maxn = nCoeffs/2;
        fourier_ns = -maxn+1:maxn;
        valuesVec = obj.apply_kron(exp(1i*obj.theta'*fourier_ns),1,fourierRect);
    end

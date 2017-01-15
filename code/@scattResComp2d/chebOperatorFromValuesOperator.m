function X = chebOperatorFromValuesOperator(obj,X)
% Input X is a map from valuesVec to valuesVec, output Y is a map from
% chebVec to chebVec.

% apply kron(U,I)*Winv on the right
X = obj.apply_Winv_on_right(obj.apply_U_kron_I_on_right(X));

% apply W*kron(Uinv,I) on the left
for col = 1:size(X,2)
    X(:,col) = obj.chebVecFromValuesVec(X(:,col));
end

function X = fourierOperatorFromValuesOperator(obj,X)
% Input X is a map from valuesVec to valuesVec, output Y is a map from
% fourierVec to fourierVec.

    X = obj.apply_U_kron_I_on_right(X);
    X = obj.apply_Uinv_kron_I(X);

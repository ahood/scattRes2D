function Y = apply_U_kron_I_on_right(obj,X)
% Returns Y = X*kron(U,I)

    Y = obj.apply_UT_kron_I(X.').';
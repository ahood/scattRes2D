function x = apply_Aorig_star(obj,x0,conjVvaluesVec)

    % need this as input to apply_kron
    xhat  = reshape(x0        , obj.Nr, []); % for applying kronecker products
    Rx = x0;
    for j = 1:size(Rx,2), Rx(:,j) = Rx(:,j)./obj.rs; end
    Rxhat = reshape(Rx, obj.Nr, []); % ditto

    x = x0;
    for j = 1:size(x,2), x(:,j) = conjVvaluesVec.*x0(:,j); end
    x = x - obj.apply_kron(obj.Dtt_star,diag(1./obj.r.^2),xhat); % Dtt term
    x = x - obj.apply_J_kron_Drr_J_star(xhat); % Drr_J term
    x = x - obj.apply_I_kron_Drr_I_star(xhat); % Drr_I term
    x = x - obj.apply_I_kron_Dr_I_star(Rxhat); % Dr_I term
    x = x - obj.apply_J_kron_Dr_J_star(Rxhat); % Dr_J term


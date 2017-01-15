function [abserr,relerr] = get_err(x_true,x_approx)
    abserr = norm(full(x_true-x_approx));
    relerr = abserr/norm(full(x_true));
end
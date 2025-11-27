function [x, a, k] = jac_iter(G, x0, tol, n)
    xr = x0;
    x = G(xr);
    a = zeros(size(x0, 1), n);
    a(:, 1) = x0;
    a(:, 2) = x;
    for k = 1:n - 1
        if norm(x - xr) < tol
            break;
        end
        xr = x;
        x = G(xr);
        a(:, k + 2) = x;
    end
    k = k + 1;
    a = a(:, 1:k + 1);
end
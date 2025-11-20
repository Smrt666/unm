function [xr, xarr, k] = iteracija(g, x0, tol, n)
    xarr = zeros(1, n);
    xarr(1) = x0;
    xarr(2) = g(x0);
    k = 1;
    while abs(xarr(k + 1) - xarr(k)) > tol && k < n
        k = k + 1;
        xarr(k + 1) = g(xarr(k));
    end
    xr = xarr(k + 1);
    xarr = resize(xarr, k + 1);
end
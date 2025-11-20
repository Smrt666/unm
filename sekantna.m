function [xr, xarr, k] = sekantna(f, x0, x1, tol, maxk)
    xarr = zeros(1, maxk);
    xarr(1) = x0;
    xarr(2) = x1;
    k = 1;
    while abs(xarr(k + 1) - xarr(k)) > tol && k < maxk
        k = k + 1;
        xarr(k + 1) = xarr(k) - f(xarr(k)) * (xarr(k) - xarr(k - 1)) / (f(xarr(k)) - f(xarr(k - 1)));
    end
    xr = xarr(k + 1);
    xarr = resize(xarr, k + 1);
end
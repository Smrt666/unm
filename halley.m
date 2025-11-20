function [x, xs, k] = halley(f, df, ddf, x0, tol, maxk)
    [x, xs, k] = iteracija(@(x) x - 2 * f(x)* df(x) / (2 * df(x)^2 - f(x) * ddf(x)), x0, tol, maxk);
end
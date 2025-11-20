function [x, xs, k] = tangentna(f, df, x0, tol, maxk)
    [x, xs, k] = iteracija(@(x) x - f(x) / df(x), x0, tol, maxk);
end
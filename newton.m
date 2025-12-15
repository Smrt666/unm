function [x, xs, k] = newton(f, j, x0, tol, n)
    x = x0 - j(x0) \ f(x0);
    xs = zeros(length(x), n);
    xs(:, 1) = x;
    k = 1;
    while norm(x - x0) > tol && k < n
        x0 = x;
        x = x - j(x) \ f(x);
        k = k + 1;
        xs(:, k) = x;
    end
    xs = resize(xs, [length(x), k]);
end
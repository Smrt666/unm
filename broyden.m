function [x, xs, k] = broyden(F, B0, x0, tol, n)
    x = x0;
    B = B0;
    xs = zeros(length(x), n);

    for k = 1:n
        s = B \ F(x);
        x_new = x - s;
        xs(:, k) = x_new;

        if norm(s) < tol
            x = x_new;
            break;
        end

        B = B + (F(x_new)*s')/(s'*s);

        x = x_new;
    end
    xs = resize(xs, [length(x), k]);
end
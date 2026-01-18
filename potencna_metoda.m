function [x, l] = potencna_metoda(A, x0, tol, n)
    x = x0;
    for r = 1:n
        x = A * x;
        x = x / norm(x);
        l = x' * A * x;
        if norm(A*x - l*x) < tol
            break
        end
    end
end
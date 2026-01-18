function [Q, R] = givens(A)
    [m, n] = size(A);
    Q = eye(m);
    for i = 1:n
        for k = i + 1:m
            r = sqrt(A(i, i)^2 + A(k, i)^2);
            c = A(i, i) / r;
            s = A(k, i) / r;
            A([i, k], i:n) = [c, s; -s, c] * A([i, k], i:n);
            Q([i, k], :) = [c, s; -s, c] * Q([i, k], :);
        end
    end
    Q = Q';
    R = A;
end
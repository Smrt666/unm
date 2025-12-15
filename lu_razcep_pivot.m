function [L, U] = lu_razcep_pivot(A)
    U = A;
    L = eye(size(A));
    n = size(A, 1);
    for i = 1:n - 1
        [~, pi] = max(maxabs((U(i:n, i:n)')));
        U([i, pi + i - 1], :) = U([pi + i - 1, i], :);
        L([i, pi + i - 1], 1:i-1) = L([pi + i - 1, i], 1:i-1);
        for j = i + 1:n
            L(j,i) = U(j,i) / U(i,i);
            for k = i:n
                U(j,k) = U(j,k) - L(j,i) * U(i,k);
            end
        end
    end
end
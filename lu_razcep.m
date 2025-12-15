function [L, U] = lu_razcep(A)
    U = A;
    L = eye(size(A));
    n = size(A, 1);
    for i = 1:n - 1
        for j = i + 1:n
            L(j,i) = U(j,i) / U(i,i);
            for k = i:n
                U(j,k) = U(j,k) - L(j,i) * U(i,k);
            end
        end
    end
end
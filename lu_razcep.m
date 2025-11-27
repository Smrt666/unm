function [L, U] = lu_razcep(A)
    L = zeros(size(A));
    U = zeros(size(A));
    for i = 1:n - 1
        for j = i + 1:n
            L(j,i) = A(j,i) / A(i,i);
            for k = i + 1:n
                U(j,k) = A(j,k) - L(j,i) * A(i,k);
            end
        end
    end
end
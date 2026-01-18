A = [4 1 -1 0;
    1 3 -1 0;
    -1 -1 5 2;
    0 0 2 4];

function [k] = obcutljivost(A, x0, tol, n)
    B = A' * A;
    [~, l] = potencna_metoda(B, x0, tol, n);
    [~, li] = potencna_metoda(inv(B), x0, tol, n);
    k = sqrt(l * li);
end

disp(obcutljivost(A, [1 1 1 1]', 0, 10));
disp(cond(A));
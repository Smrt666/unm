A = [4 1 -1 0;
    1 3 -1 0;
    -1 -1 5 2;
    0 0 2 4];

[x, l] = potencna_metoda(A, [1;1;1;1], 0, 10);
disp(l)
eig(A)

function [x, l] = potencna_metoda2(f, x0, tol, n)
    x = x0;
    for r = 1:n
        x = f(x);
        x = x / norm(x);
        l = x' * f(x);
        if norm(f(x) - l*x) < tol
            break
        end
    end
end

[x, l] = potencna_metoda2(@(x) A \ x, [1;1;1;1], 0, 10);
disp(1/l);

A = rosser;
B = A;
B(1,8) = -29;
C = A(:,8:-1:1);

function [M] = qr_iteracija(A, n)
    M = A;
    for i = 1:n
        [Q, R] = qr(M);
        M = R * Q;
    end
end


diag(qr_iteracija(A, 100));



diag(qr_iteracija(hess(A), 100));


function [lr] = qr_iteracija2(A, tol, n)
    M = A;
    if isscalar(A)
        lr = A;
        return
    end
    for i = 1:n
        if abs(M(end, end-1)) < (abs(M(end-1, end-1)) + abs(M(end, end))) * tol
            break
        end
        [Q, R] = qr(M);
        M = R * Q;
    end
    lr = [M(end, end); qr_iteracija2(M(1:end-1, 1:end-1), tol, n)];
end

disp(qr_iteracija2(A, 1e-10, 1e6));

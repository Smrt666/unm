n = 1000;
A = spdiags([1 -2 1],-1:1,n,n);

% disp(A);
tic
a2 = normest(A);
toc
disp(a2);
disp(a2 - norm(full(A), 2));
tic
disp(sqrt(max(eig(A' * A))));
toc
tic
disp(max(svd(A)));
toc

tic
a1 = norm(A, 1);
ainf = norm(A, inf);
af = norm(A, "fro");
aninf = norm(reshape(A, [], 1), inf);

sqn = sqrt(n);
lb = max([af / sqn, ainf / sqn, a1 / sqn, aninf]);
ub = min([af, sqn * ainf, sqn * a1, n * aninf, sqrt(a1 * ainf)]);
toc

fprintf('Lower bound: %.4f\nUpper bound: %.4f\n', lb, ub);
disp(a2);
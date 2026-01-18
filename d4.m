x = (1/5 * (0:10))';
f = @(x) x .* sin(3.*x);
y = f(x);

degree = 3;
A = zeros(length(x), degree + 1);
for i = 0:degree
    A(:, i + 1) = x .^ i;
end

a = A \ y;
disp(a);

degree = 3;
x2 = x([1:5, 7:11]);
y2 = f(x2) - f(1);
K = [x2 - 1, x2.^2 - 1, x2.^3 - 1];

disp(K);
disp(y2);

k = K \ y2;
k = [f(1) - sum(k); k];
disp(k);
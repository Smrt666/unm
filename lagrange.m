function [p] = lagrangeev_polinom(x, i)
    n = length(x);
    p = 1;
    i = i + 1;
    for k = 1:n
        if k ~= i
            p = conv(1 / (x(i) - x(k)) * [1, -x(k)], p);
        end
    end
end

x = [0, 1/2, 3/2, 3, 5];
for i = 0:4
    disp(polyval(lagrangeev_polinom(x, i), 2));
end

f = @(x) x .* sin(pi * x);
p = [0, 0, 0, 0, 0];
for i = 1:5
    p = p + f(x(i)) * lagrangeev_polinom(x, i - 1);
end

X = linspace(0, 5);
Y = polyval(p, X);
plot(X, Y);
xlabel('X');
ylabel('Interpolated Values');
title('Lagrange Interpolation');

hold on;
Yf = f(X);
plot(X, Yf, '--');
legend('Lagrange Interpolation', 'Original Function');
hold off;
f = @(x) x + 4 - exp(x^2);
df = @(x) 1 - 2*x*exp(x^2);
ddf = @(x) 2 * exp(x^2) + 4 * x^2 * exp(x^2);

tol = 1e-10;
n = 100;

function [kr] = mink(xs, x0)
    for i = 1:length(xs)
        if abs(xs(i) - x0) < 1e-10
            kr = i;
            break
        end
    end
    kr = length(xs);
end

x0s = 1:0.2:10;
tank = zeros(1, length(x0s));
hallk = zeros(1, length(x0s));
sekk = zeros(1, length(x0s));
for i = 1:length(x0s)
    x0 = x0s(i);
    zr = fzero(f, x0);
    [~, xs, ~] = tangentna(f, df, x0, tol, n);
    tank(i) = mink(xs, zr);
    [~, xs, ~] = halley(f, df, ddf, x0, tol, n);
    hallk(i) = mink(xs, zr);
    [~, xs, ~] = sekantna(f, x0, x0 + 0.1, tol, n);
    sekk(i) = mink(xs, zr);
end

disp(x0s);
disp(tank);

hold on
scatter(x0s, tank);
scatter(x0s, hallk);
scatter(x0s, sekk);
hold off
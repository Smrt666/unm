px = [-17, 52, -42, 8, 0];
py = [-19.5, 64, -60, 16, 0];

t = linspace(0, 1);

x = polyval(px, t);
y = polyval(py, t);
plot(x, y);
hold on
inside = (x .* x + y .* y > 0.25) & (x .* x + y .* y < 1);

x(~inside) = NaN;
y(~inside) = NaN;
plot(x(inside), y(inside));
plot(cos(2 * pi * t), sin(2 * pi * t));

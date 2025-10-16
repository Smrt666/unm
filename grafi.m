figure
x1 = linspace(1, 3);
y1 = sin(x1) .* exp(sqrt(x1));
plot(x1, y1);

figure
t = linspace(0, 2*pi);
x2 = sin(t);
y2 = cos(t);
plot(x2, y2);

figure
[x3, y3] = meshgrid(linspace(0, 1), linspace(0,1));
z3 = (x3.^2 + y3.^2) ./ (1 + x3 + y3);
surf(x3, y3, z3);


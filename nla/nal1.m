xy = [
    41.3477   41.1652   40.9966   40.8009   40.6111   40.4934   40.2494   40.0246   39.8930   39.6418   39.4306   39.1946   39.0827   38.8907   38.7291;
    -17.8104  -17.7921  -17.7105  -17.7045  -17.7041  -17.7034  -17.7548  -17.8773  -17.9660  -18.1852  -18.2932  -18.4563  -18.6589  -18.8606  -19.1617
];

x = xy(1, :);
y = xy(2, :);

disp(x);
A = [x.*y; y.*y; x; y; ones(size(x))]';
target = [-x.*x]';

bcde = A \ target;
abcde = [1; bcde];

disp(bcde(1))

fh = @(xx,yy) (abcde' * [xx.*xx; yy.*xx; yy.*yy; xx; yy; ones(size(xx))]);

axis equal
grid on
fimplicit(fh);

% Sestavimo matriko A za TLS
A = [x(:).*y(:), y(:).^2, x(:), y(:), ones(length(x),1), x(:).^2];

% TLS rešitev preko SVD
[~, ~, V] = svd(A, 0);
p = V(:, end);   % lastni vektor za najmanjšo singularno vrednost

% Parametri (normaliziramo tako, da je a = 1)
a = p(6);
b = p(1)/a;
c = p(2)/a;
d = p(3)/a;
e = p(4)/a;
f = p(5)/a;

% =========================
% 2) FUNKCIJA ELIPSE y = p(x)
% (rešimo kvadratno enačbo po y)
% =========================
ellipse_y = @(xq) ...
    (- (b*xq + e) + sqrt((b*xq + e).^2 - 4*c*(xq.^2 + d*xq + f))) ./ (2*c);

% =========================
% 3) RAZDALJA DO TOČKE T
% =========================
T = [38.8907; -18.8606];

dist = @(xq) norm([xq; ellipse_y(xq)] - T);

% =========================
% 4) MINIMIZACIJA (fminbnd)
% =========================
xmin = min(x);
xmax = max(x);

[x_opt, d_min] = fminbnd(dist, xmin, xmax);

y_opt = ellipse_y(x_opt);

fprintf('Najbližja točka na elipsi: (%.4f, %.4f)\n', x_opt, y_opt);
fprintf('Najmanjša razdalja: %.15f\n', d_min);

% =========================
% 5) RISANJE
% =========================
fimplicit(@(x,y) x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f, ...
          [xmin xmax min(y) max(y)], 'LineWidth', 1.5);
hold on;
plot(T(1), T(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(x_opt, y_opt, 'bx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Elipsa (TLS)', 'Točka T', 'Najbližja točka');
grid on;
axis equal;

% TLS
% C = [A b];
% [U,S,V] = svd(C,0);
% v = V(:,end);
% 
% x_tls = -v(1:end-1)/v(end);
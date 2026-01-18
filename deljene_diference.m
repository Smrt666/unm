function [n, d] = deljene_diferencee(x, y)
% n: newtonovi keficienti (baza: 1, (x-x0), (x-x0)(x-x1), ...)
% d: shema diferenc
    d = zeros(length(x));
    d(:, 1) = y;
    for j = 2:length(x)
        d(1:end-j+1, j) = (d(2:end - j + 2, j - 1) - d(1:end - j + 1, j - 1)) ./ (x(j:end) - x(1:end - j + 1))';
    end
    n = d(1, :);
end

function [p] = newtonov_bazni(x, n)
    p = 1;
    for k = 1:n
        p = conv(p, [1, -x(k)]);
    end
end

function [p] = newtonov_zapis(x, y)
    np = deljene_diferencee(x, y);
    n = length(x);
    p = 0;
    for i = 0:n - 1
        p = p + [zeros(1, n - i), np(i + 1) * newtonov_bazni(x, i)];
    end
end

f = @(x) 40 ./ (x + 1);
x = [0, 1/3, 2/3, 1];

n = deljene_diferencee(x, f(x));
disp(n);

p = newtonov_zapis(x, f(x));
disp(p);
disp(polyfit(x, f(x), length(x) - 1));

x = linspace(0, 1, 1001);
disp(max(abs(polyval(p, x) - f(x))));

function [c, D] = deljene_diference_z_odvodom(x, fx, dfx)
    n = length(x);
    m = 2 * n;

    z = zeros(1, m);
    D = zeros(m, m);

    % Podvojimo točke
    for i = 1:n
        z(2 * i - 1) = x(i);
        z(2 * i) = x(i);

        D(2 * i - 1, 1) = fx(i);
        D(2 * i, 1) = fx(i);
    end

    % Začetni stolpec deljenih diferenc
    for i = 1:m-1
        if mod(i,2) == 1
            D(i, 2) = dfx((i + 1) / 2);
        else
            D(i, 2) = (D(i + 1, 1) - D(i, 1)) / (z(i + 1) - z(i));
        end
    end

    % Ostale deljene diference
    for j = 3:m
        for i = 1:m-j+1
            D(i, j) = (D(i + 1, j - 1) - D(i, j - 1)) / (z(i + j - 1) - z(i));
        end
    end
    c = D(1, :);
end

function [p] = newtonov_zapis2(x, np)
    n = length(x);
    p = 0;
    for i = 0:n - 1
        p = p + [zeros(1, n - i), np(i + 1) * newtonov_bazni(x, i)];
    end
end

% a: Seznam deljenih diferenc
% x: Seznam premaknjenih potenc
% t: Točka, v kateri izračunamo polinom
function y = newtonov_polinom(a, x, t)
    n = length(a);
    y = a(1) * ones(size(t));
    prod = ones(size(t));

    for k = 2:n
        prod = prod .* (t - x(k-1));
        y = y + a(k) * prod;
    end
end

x = [0, 1];
y = [40, 20];
dy = [-40, -10];

c = deljene_diference_z_odvodom(x, y, dy);
cp = newtonov_zapis2([0, 0, 1, 1], c);
x = linspace(0, 1, 1001);
disp(max(abs(polyval(cp, x) - f(x))));
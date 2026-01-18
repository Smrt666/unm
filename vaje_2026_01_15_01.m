% Izračun deljenih diferenc in vrednosti Newtonovega polinoma

function D = deljene_diference_brez_odvoda(x, fx)
    n = length(x);

    D = zeros(n);
    D(:, 1) = fx;

    for j = 2:n
        for i = 1:n-j+1
            D(i, j) = (D(i + 1, j - 1) - D(i, j - 1)) / (x(i + j - 1) - x(i));
        end
    end
end

function D = deljene_diference_z_odvodom(x, fx, dfx)
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

% Uporaba

f = @(x) 40 ./ (x + 1);
df = @(x) -40 ./ (x + 1).^2;

% Brez odvoda

x0 = [0, 1/3, 2/3, 1];
fx0 = f(x0);

D0 = deljene_diference_brez_odvoda(x0, fx0); % Deljene diference brez odvoda
a0 = D0(1, :); % Newtonovi koeficienti brez odvoda

xt0 = (0:1000) / 1000; % Vrednosti za izračun newtonovega polinoma
pt0 = newtonov_polinom(a0, x0, xt0); % Vrednosti newtonovega polinoma

err0 = max(abs(f(xt0) - pt0));

% Z odvodom

x1 = [0, 1];
fx1 = f(x1);
dfx1 = df(x1);

D1 = deljene_diference_z_odvodom(x1, fx1, dfx1); % Deljene diference z odvodom
a1 = D1(1, :); % Newtonovi koeficienti z odvodom
disp(a1);

z1 = repelem(x1, 2); % Podvojimo premaknjene potence

xt1 = (0:1000) / 1000; % Vrednosti za izračun newtonovega polinoma
pt1 = newtonov_polinom(a1, z1, xt1); % Vrednosti newtonovega polinoma

err1 = max(abs(f(xt1) - pt1));

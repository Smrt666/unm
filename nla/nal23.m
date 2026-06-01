% 0. Priprava mreže in generiranje podatkov
N = 100;
t = linspace(-5, 5, N);
[X2, Y2] = ndgrid(t, t); % Za matriko M1 (2D)
[X3, Y3, Z3] = ndgrid(t, t, t); % Za tenzor M2 (3D)

% f1(x,y)
M1 = 1.4000 * cos(X2).^2 + 1.1000 * sin(Y2).^2 + X2.*Y2 + sin(X2.*Y2 / 5.0);

% f2(x,y,z)
M2 = 1.4000 * cos(X3).^2 + 1.1000 * sin(Y3).^2 + X3.*Y3 + sin(X3.*Y3 / 5.0) + 0.10 * Z3.^2;

% 1. DEL: Rang m in absolutna razlika
[U, S, V] = svd(M1);
max_err = Inf;
m = 0;
M_m = zeros(N, N);

while max_err >= 1e-3
    m = m + 1;
    % Dodamo m-ti komponento k aproksimaciji
    M_m = M_m + S(m,m) * U(:, m) * V(:, m)';
    % Izračun maksimalne absolutne razlike
    max_err = max(abs(M1(:) - M_m(:)));
end

fprintf('--- 1. DEL ---\n');
fprintf('Potreben rang m: %d\nAbsolutna razlika: ', m);
disp(max_err)


% 2. DEL: Tenzorske operacije
% P = M2 x_2 M1

P = ttm(tensor(M2), M1, 2);

pr = double(tenmat(P, 3));

% Izračun druge norme 5000. stolpca matrike P_tau
stolpec_5000 = pr(:, 5000);
druga_norma = norm(stolpec_5000, 2);

fprintf('--- 2. DEL ---\n');
fprintf('Druga norma 5000. stolpca matrike P_r: %.12f\n', druga_norma);

% cp_als(x, r) x na r komponent (obstaja tut tucker_als)
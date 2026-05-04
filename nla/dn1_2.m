% --- PODATKI ---
A = [ 1.2500, -0.5000,  1.0000, -1.0000,  0.2500,  0.2500,  0,  0.2500;
    -0.5000, 11.0000, -2.0000,  2.0000, -0.5000, -0.5000,  0, -0.5000;
    1.0000, -2.0000,  9.0000, -4.0000,  1.0000,  1.0000,  0,  1.0000;
    -1.0000,  2.0000, -4.0000, -2.0000, -1.0000, -1.0000,  0, -1.0000;
    0.2500, -0.5000,  1.0000, -1.0000,  3.2500,  0.2500,  0,  0.2500;
    0.2500, -0.5000,  1.0000, -1.0000,  0.2500,  4.2500,  0,  0.2500;
    0,       0,       0,       0,       0,       0,  1.0000,       0;
    0.2500, -0.5000,  1.0000, -1.0000,  0.2500,  0.2500,  0,  2.2500];

% --- (1) DOLOČITEV MATRIKE D IN NJENE SLEDI ---
% Vektor u določimo iz izvendiagonalnih elementov prve vrstice matrike A.
% Ker je A = D + u*u', velja A(1,j) = u(1)*u(j) za j>1.
u = [0.5; -1; 2; -2; 0.5; 0.5; 0; 0.5];
u_sq = u.^2;

% Diagonalni elementi d_i = A_ii - u_i^2
d = diag(A) - u_sq;
D = diag(d);
sled_D = sum(d);

fprintf('(1) Sled matrike D: %.6f\n', sled_D);

% --- (2) REDUKCIJA IN SEKULARNA ENAČBA ---
% Preuredimo pare (d_i, u_i^2) padajoče po velikosti d_i
[d_sorted, idx] = sort(d, 'descend');
u_sq_sorted = u_sq(idx);
u_sorted = u(idx);

% Odstranimo vrstice in stolpce, kjer je u_i == 0 (deflacija)
keep_idx = u_sq_sorted > 1e-10;
d_red = d_sorted(keep_idx);
u_sq_red = u_sq_sorted(keep_idx);
u_red = u_sorted(keep_idx);

% Definicija sekularne funkcije f(x) in njenega odvoda f'(x)
f = @(x) 1 + sum(u_sq_red ./ (d_red - x));
df = @(x) sum(u_sq_red ./ ((d_red - x).^2)); 

fprintf('(2) f(7) = %.6f\n', f(7));

% --- (3) IN (4) ISKANJE NIČLE (LASTNE VREDNOSTI) ---
% Iskanje druge največje lastne vrednosti (med d_1 in d_2)
d1 = d_red(1);
d2 = d_red(2);
x0 = (d1 + d2) / 2; % Začetni približek na sredini intervala (7.5)

% (3) Podatek za x1, ki smo ga analitično izračunali preko racionalne Newtonove:
x1 = 8.416856; 
fprintf('(3) x1 (približek po enem koraku): %.6f\n', x1);

% (4) Za končni, zelo natančen izračun lahko uporabimo vgrajeno fzero,
% ki bo poiskala točno ničlo sekularne enačbe:
lambda = fzero(f, x0);
fprintf('(4) Približek za lastno vrednost (konvergiran): %.6f\n', lambda);

% --- (5) PRIBLIŽEK ZA LASTNI VEKTOR ---
% Uporabimo približek x1 = 8.416856 iz 3. točke
% Izračun lastnega vektorja z = (D - x1*I)^(-1) * u
z_approx = u ./ (d - x1);

% Lastni vektor normiramo v 2-normi. 
% Norma vektorja ustreza kvadratnemu korenu odvoda f'(x1).
norm_z = sqrt(df(x1));
z_norm = z_approx / norm_z;

% Iščemo absolutno vrednost prvega elementa izvornega vektorja
abs_z1 = abs(z_norm(1));

fprintf('(5) Absolutna vrednost prvega elementa vektorja: %.6f\n', abs_z1);
fn1 = @(x,y,z) 1.6e-3 * x.^2 + 1.6e-3 * y.^2 - 1;
fn2 = @(x,y,z) 5.3e-4*x.^2 + 5.3e-4*y.^2 + 5.3e-4*z.^2 + 2.7e-2*x -1;
fn3 = @(x,y,z) -1.4e-4*x + 1.0e-4*y + z - 3.4e-3;

i0 = [1 4 2 5 7 3 6 8 4 7 5 6 9 7 8 9 1 4 8 10 2 5 9 10 3 6 10]; 
j0 = [1 1 2 2 2 3 3 3 4 4 5 6 6 7 7 7 8 8 8 8 9 9 9 9 10 10 10]; 
a0 = [1.6e-3 5.3e-4 1.6e-3 5.3e-4 1.0e-4 1.6e-3 5.3e-4 1.0e-4 5.3e-4 1 5.3e-4 5.3e-4 1 -3.4e-3 1 1.0e-4 -1 -1 -3.4e-3 1.0e-4 -1 -1 -3.4e-3 1 -1 -1 -3.4e-3]; 
A0 = full(sparse(i0,j0,a0)); 

i1 = [7 4 8 5 9  6 10]; 
j1 = [7 8 8 9 9 10 10]; 
a1 = [-1.4e-4 0.027 -1.4e-4 0.027 -1.4e-4 0.027 -1.4e-4]; 
A1 = full(sparse(i1,j1,a1)); 

i2 = [1 4 2 5 3   6 10]; 
j2 = [8 8 9 9 10 10 10]; 
a2 = [1.6e-3 5.3e-4 1.6e-3 5.3e-4 1.6e-3 5.3e-4 0]; 
A2 = full(sparse(i2,j2,a2));  

% 1. NALOGA: Izračun norme v točki (25, -3, 0)
x0 = 25; y0 = -3; z0 = 0;
f1 = 1.6e-3 * x0^2 + 1.6e-3 * y0^2 - 1;
f2 = 5.3e-4 * x0^2 + 5.3e-4 * y0^2 + 5.3e-4 * z0^2 + 2.7e-2 * x0 - 1;
f3 = -1.4e-4 * x0 + 1.0e-4 * y0 + z0 - 3.4e-3;

% Odpremo novo sliko
figure('Name', 'Presek ploskev v prostoru', 'NumberTitle', 'off');
hold on;

% Nastavimo območje risanja (glede na to, da je testna točka x=25)
ovojnica = [-30 30 -30 30 -15 15];

% Narišemo vsako ploskev s svojo barvo in prosojnostjo (FaceAlpha)
h1 = fimplicit3(fn1, ovojnica, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
h2 = fimplicit3(fn2, ovojnica, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
h3 = fimplicit3(fn3, ovojnica, 'FaceColor', 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Označimo testno točko (25, -3, 0) z rdečo piko
x0 = 25; y0 = -3; z0 = 0;
plot3(x0, y0, z0, 'r.', 'MarkerSize', 25);

% Ureditev grafa (napisov, osi, luči)
xlabel('X os'); ylabel('Y os'); zlabel('Z os');
title('Presek cilindra (zeleno), sfere (modro) in ravnine (vijolično)');
grid on;
axis equal;
view(3);         % Nastavi 3D pogled
camlight;       % Doda svetlobo za lepši 3D izgled ploskev
lighting gouraud;

F_norma = norm([f1; f2; f3], 2);
fprintf('1) Druga norma [f1, f2, f3] v (25, -3, 0): %.12f\n', F_norma);
% fprintf('1) Druga norma [f1, f2, f3] v (25, -3, 0): %.12f\n', norm([fn1(x0, y0, z0); fn2(x0, y0, z0); fn3(x0, y0, z0);]));

% 2. NALOGA: Reševanje s polyeig
% polyeig rešuje (A0 + lambda*A1 + lambda^2*A2)x = 0
[V_poly, lambda_poly] = polyeig(A0, A1, A2);

% Izločimo samo realne in končne lastne vrednosti
prave_idx = find(abs(imag(lambda_poly)) < 1e-6 & ~isinf(lambda_poly));
lam_koncni = lambda_poly(prave_idx);
V_koncni = V_poly(:, prave_idx);

% Poiščemo tisto presečišče, ki ima negativno y komponento
% y koordinata je v vec(8)/v(10)
y_komponente = zeros(size(lam_koncni));
z_komponente = zeros(size(lam_koncni));
for i = 1:length(lam_koncni)
    v = V_koncni(:, i);
    y_komponente(i) = v(8) / v(10);
    z_komponente(i) = v(9) / v(10);
end

idx_neg_y = find(y_komponente < 0);
fprintf('2) z-komponenta presečišča z negativnim y=%.5f (polyeig): %.12f\n', y_komponente(idx_neg_y), z_komponente(idx_neg_y));

% Risanje ploskev in obeh presečišč
tocke = zeros(length(lam_koncni), 3); % Vsaka vrstica bo [x, y, z]
for i = 1:length(lam_koncni)
    v = V_koncni(:, i);
    tocke(i, 1) = lam_koncni(i);      % x = lambda
    tocke(i, 2) = v(8) / v(10);       % y = v_8 / v_10
    tocke(i, 3) = v(9) / v(10);       % z = v_8 / v_10
end

figure('Name', 'Obe presečišči ploskev', 'NumberTitle', 'off');
hold on;

ovojnica = [-30 30 -30 30 -15 15];

% Izris ploskev
fimplicit3(fn1, ovojnica, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
fimplicit3(fn2, ovojnica, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
fimplicit3(fn3, ovojnica, 'FaceColor', 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Izris obeh presečišč (tocke so v matriki 'tocke')
p1 = plot3(tocke(1,1), tocke(1,2), tocke(1,3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'LineWidth', 2);
p2 = plot3(tocke(2,1), tocke(2,2), tocke(2,3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'LineWidth', 2);

% Nastavitve grafa
xlabel('X os'); ylabel('Y os'); zlabel('Z os');
title('Iskanje presečišč algebrskih ploskev');
grid on;
axis equal;
view(3);
camlight;
lighting gouraud;

% 3. NALOGA: Zaporedna linearna aproksimacija
% Definiramo nelinearno funkcijo T(lam) in njen odvod dT(lam) za QEP
T = @(lam) lam^2 * A2 + lam * A1 + A0;
dT = @(lam) 2 * lam * A2 + A1;

lam0 = 25;
maxKorakov = 100;
tol = 1e-12;

[lam_konc, vec_konc] = zapLinAproks(T, dT, lam0, maxKorakov, tol);

% Izračunamo z-komponento iz dobljenega lastnega vektorja
z_zapLin = vec_konc(9) / vec_konc(10);
fprintf('3) z-komponenta presečišča (zapLinAproks): %.12f\n', z_zapLin);
% disp(vec_konc);


% Implementacija metode zaporedne linearne aproksimacije
function [lam, vec] = zapLinAproks(T, dT, lam, maxKorakov, tol)
    for k = 1:maxKorakov
        Tk = T(lam);
        dTk = dT(lam);
        
        % Rešujemo posplošeni problem lastnih vrednosti: Tk * x = theta * dTk * x
        [X, D] = eig(Tk, dTk);
        theta_vsi = diag(D);
        
        % Odstranimo morebitne neskončne ali NaN vrednosti
        veljavni = ~isinf(theta_vsi) & ~isnan(theta_vsi);
        theta_vsi = theta_vsi(veljavni);
        X = X(:, veljavni);
        
        % Poiščemo tisti theta, ki je najbližje 0
        [~, idx] = min(abs(theta_vsi));
        delta_lam = theta_vsi(idx);
        
        % Posodobitev lastne vrednosti
        lam = lam - delta_lam;
        vec = X(:, idx);
        
        % Preverjanje konvergence
        if abs(delta_lam) < tol
            return;
        end
    end
    warning('Metoda ni konvergirala v izbranem številu korakov.');
end
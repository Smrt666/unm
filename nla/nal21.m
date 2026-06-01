function [R2, Q] = zamenjaj(R, i, j)
    % ZAMENJAJ Zamenja diagonalca R(i,i) in R(j,j) v realni zgornje trikotni matriki.
    % Vsak korak sosednje zamenjave uporabi par Givensovih rotacij, ki ohranita
    % zgornje trikotno strukturo in zadoščata pogoju R2 = Q * R * Q'.
    
    n = size(R, 1);
    Q = eye(n);
    R2 = R;
    
    if i == j, return; end
    
    % Vedno uredimo indekse tako, da je i < j
    if i > j
        tmp = i; i = j; j = tmp;
    end
    
    % 1. i-ti element potiskamo DESNO do mesta j
    for k = i : (j - 1)
        [R2, Gk] = zamenjaj_sosednja(R2, k);
        Q = Gk * Q;
    end
    
    % 2. Prvotni j-ti element potiskamo LEVO nazaj do mesta i
    for k = (j - 2) : -1 : i
        [R2, Gk] = zamenjaj_sosednja(R2, k);
        Q = Gk * Q;
    end
end

function [R_new, Q] = zamenjaj_sosednja(R, i)
    %   R - zgornja trikotna matrika (n x n)
    %   i - indeks, kjer želimo zamenjati i-ti in (i+1)-ti diagonalni element
    %
    % Izhoda:
    %   R_new - nova zgornja trikotna matrika z zamenjanima elementoma
    %   Q     - ortogonalna transformacijska matrika (Q * R * Q' = R_new)

    n = size(R, 1);
    
    % Začnemo z identično matriko
    Q = eye(n);
    
    % Izvlečemo vrednosti iz 2x2 bloka na diagonali
    a = R(i, i);
    b = R(i, i+1);
    c = R(i+1, i+1);
    
    x = [b; c - a];
    
    % planerot vrne 2x2 Givensovo rotacijo G, takšno da je G*x = [norm(x); 0]
    [G, ~] = planerot(x);
    
    % Vstavimo 2x2 rotacijo v ustrezno mesto n x n matrike Q
    Q(i:i+1, i:i+1) = G;
    
    % Izračunamo novo matriko (podobnostna transformacija)
    R_new = Q * R * Q';
    
    % Zaradi numeričnih napak pri računanju s plavajočo vejico eksaktno
    % postavimo element pod diagonalo na 0, da ohranimo trikotno strukturo.
    R_new(i+1, i) = 0;
end

A = magic(6);
[U,R] = schur(A);

[R2, Q] = zamenjaj(R, 3, 4);

disp(norm(R - R2, 2));
disp(norm(Q, 'fro'));

[R3, Q2] = zamenjaj(R, 1, 6);

disp(norm(R - R3, 2));
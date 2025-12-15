function [L, U, P, Q] = lu_kompletni(A)     %Določimo matrike L, U, P, Q, da je P*A*Q = L*U
[n,~] = size(A);    %n = velikost matrike
P = eye(n);         %inicializiramo P
Q = eye(n);         %inicializiramo Q
for j = 1 : n-1
    AA = abs(A(j:n, j:n))';     %Podmatrika |A|^T od j-tega do n-tega stolpca in od j-te do n-te vrstice
    [arr, b] = max(AA);         %arr...max elementi vsake vrstice v |A|, b...indeksi, pripadajoči max elementom
    [x, i] = max(arr);          %x...max element v AA, i pripadajoči indeks
    max_vrstica = i + j - 1;    %vrstica, ki vsebuje element x
    max_stolpec = b(i) + j - 1; %stolpec, ki vsebuje element x
    A([j max_vrstica], :) = A([max_vrstica j], :);  %zamenjava vrstic v A
    P([j max_vrstica], :) = P([max_vrstica j], :);  %zamenjava vrstic v P
    A(:, [j max_stolpec]) = A(:, [max_stolpec j]);  %zamenjava stolpcev v A. Tako max element pripelješ na ustrezno mesto.
    Q(:, [j max_stolpec]) = Q(:, [max_stolpec j]);  %zamenjava stolpcev v Q
    for i = j+1 : n     %Standardni postopek za LU razcep
        A(i,j) = A(i,j)/A(j,j);
        for k = j+1 : n
            A(i, k) = A(i, k) - A(i, j)*A(j, k);
        end
    end
end
L = tril(A, -1) + eye(n);   %Za L vzameš elemente A-ja pod diagonalo, na diagonalo daš enke
U = triu(A);        %Za U vzameš zgornji trikotnik A-ja 
end
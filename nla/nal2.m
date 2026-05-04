A = [1.2500   -0.5000    1.0000   -1.0000    0.2500    0.2500         0    0.2500;
    -0.5000   11.0000   -2.0000    2.0000   -0.5000   -0.5000         0   -0.5000;
    1.0000   -2.0000    9.0000   -4.0000    1.0000    1.0000         0    1.0000;
    -1.0000    2.0000   -4.0000   -2.0000   -1.0000   -1.0000         0   -1.0000;
    0.2500   -0.5000    1.0000   -1.0000    3.2500    0.2500         0    0.2500;
    0.2500   -0.5000    1.0000   -1.0000    0.2500    4.2500         0    0.2500;
    0         0         0         0         0         0    1.0000         0;
    0.2500   -0.5000    1.0000   -1.0000    0.2500    0.2500         0    2.2500
];

% u = [1, -2, 4, -4, 1, 1, 0, 1];
% disp(A - u' * u / 4);

function [D, p, u] = decompose_A(A)
    % Nisem in ne mislm brat tega
    n = size(A, 1);
    p_u_kvadrat = zeros(n, 1);
    
    % 1. Poiščemo indekse, kjer u_i ni 0 (tisti, ki imajo ne-ničelne izven-diagonale)
    A_off = A - diag(diag(A));
    ind_non_zero = find(any(A_off ~= 0, 2));
    
    if length(ind_non_zero) < 3
        warning('Premalo ne-ničelnih elementov v u za enolično rekonstrukcijo.');
    end
    
    % 2. Izračunamo p*u_i^2 samo za ne-ničelne vrstice
    for i = 1:n
        if ismember(i, ind_non_zero)
            % Najdemo dva indeksa j, k, kjer so A(i,j), A(i,k) in A(j,k) != 0
            % To storimo tako, da izberemo stolpce z največjo energijo
            ostali = setdiff(ind_non_zero, i);
            najdeno = false;
            for idx_j = 1:length(ostali)
                for idx_k = (idx_j + 1):length(ostali)
                    j = ostali(idx_j);
                    k = ostali(idx_k);
                    if A(j,k) ~= 0
                        p_u_kvadrat(i) = (A(i,j) * A(i,k)) / A(j,k);
                        najdeno = true;
                        break;
                    end
                end
                if najdeno, break; end
            end
        else
            % Če je vrstica prazna, je u_i = 0
            p_u_kvadrat(i) = 0;
        end
    end
    
    % 3. Rekonstrukcija D
    D = diag(diag(A) - p_u_kvadrat);
    
    % 4. Rekonstrukcija p in u (normalizirano na ||u|| = 1)
    p = sum(p_u_kvadrat);
    if p == 0
        u = zeros(n, 1);
    else
        u_abs = sqrt(p_u_kvadrat / p);
        u = zeros(n, 1);
        
        % Določimo znake glede na prvi ne-ničelni element
        prvi = ind_non_zero(1);
        u(prvi) = u_abs(prvi);
        for i = 1:n
            if i == prvi
                continue;
            elseif ismember(i, ind_non_zero)
                % Znak določi A(prvi, i) = p * u(prvi) * u(i)
                if sign(A(prvi, i)) ~= sign(p * u(prvi))
                    u(i) = -u_abs(i);
                else
                    u(i) = u_abs(i);
                end
            end
        end
    end
end

function [D, P, U] = rekonstruiraj_svd_iter(A, rank_k)
    % Tole je nekej kar je AI skuhow pa dela.
    % Dela tut nekak če bi hotu obliko A = D + p u u' + r v v' + ...,
    % (prvih rank_k stolpcev U so u, v, ...; p, r sta v P)
    % ampak konvergenca je za rank_k > 1 zelo slaba
    n = size(A, 1);
    
    % Inicializacija: diagonalo za začetek postavimo na povprečje izven-diagonalnih
    % To pomaga pri konvergenci, če je diagonala v D zelo velika.
    L = A;
    off_diag_avg = mean(abs(A(~eye(n))));
    for j = 1:n
        L(j,j) = off_diag_avg; 
    end
    
    for iter = 1:100
        % 1. Najboljši nizkorangovni približek (rank-k)
        [U, S, V] = svds(L, rank_k);
        L_rank = U * S * V';
        
        % 2. Posodobimo diagonalo:
        % Diagonala L mora postati enaka diagonali L_rank, 
        % izven-diagonalni elementi pa morajo ostati enaki tistim v A.
        for j = 1:n
            L(j,j) = L_rank(j,j);
        end
    end
    
    % Končni rezultati
    D_diag = diag(A) - diag(L_rank);
    D = diag(D_diag);
    U = U; % Matrika z ortonormiranimi stolpci
    P = diag(S); % Vektor uteži (p_i)
end

[D, p, u] = decompose_A(A);
d = diag(D);
disp(u');

disp(sum(d));

[d, idx] = sort(d, 'descend');
u = u(idx);
disp(d');

f = @(x) 1 + sum(p * u .* u ./ (d - x));
disp(f(7));

function x = racionalni_newton(d, u, ro, k, i, x0)
    % če je ro < 0 je največja lastna vrednost med d1 in d2
    % če je ro > 0 je največja lastna vrednost med +inf in d1
    % podobno za majmanjšo

    % d: vektor polov (urejen padajoče: d1 > d2 > ... > dn)
    % u: vektor uteži
    % ro: skalar (p)
    % k: število korakov
    % i: indeks intervala (med d(i+1) in d(i))
    % x0: (opcijsko) začetni približek
    
    % Preveri začetni približek
    if nargin < 6 || isempty(x0)
        x = (d(i) + d(i+1)) / 2;
    else
        x = x0;
    end
    
    n = length(d);
    
    for step = 1:k
        % Razdelimo sumande na tiste desno (psi) in tiste levo (phi) od intervala
        % j = 1:i so poli >= d(i) (desno od x)
        % j = i+1:n so poli <= d(i+1) (levo od x)
        
        idx_desno = 1:i;
        idx_levo = (i+1):n;
        
        % Izračun delnih vsot in njihovih odvodov
        psi_val = sum((ro * u(idx_desno).^2) ./ (d(idx_desno) - x));
        phi_val = sum((ro * u(idx_levo).^2) ./ (d(idx_levo) - x));
        
        psi_der = sum((ro * u(idx_desno).^2) ./ (d(idx_desno) - x).^2);
        phi_der = sum((ro * u(idx_levo).^2) ./ (d(idx_levo) - x).^2);
        
        % Vrednost funkcije f(x)
        fx = 1 + psi_val + phi_val;
        
        % Koeficienti racionalne aproksimacije
        % Aproksimiramo: f(y) \approx 1 + a + b/(d(i)-y) + c + d_const/(d(i+1)-y)
        % Opomba: d_const pišem kot 'd_val', da ni zmede z vektorjem d
        
        delta_i = d(i) - x;        % razdalja do desnega pola
        delta_i1 = d(i+1) - x;     % razdalja do levega pola (negativna)
        
        % b in d_val določimo iz odvodov najbližjih polov
        b = psi_der * delta_i^2;
        d_val = phi_der * delta_i1^2;
        
        % a in c določimo tako, da se ujemata vrednosti psi in phi
        a = psi_val - b/delta_i;
        c = phi_val - d_val/delta_i1;
        
        % Rešujemo: (1 + a + c) + b/(d(i)-y) + d_val/(d(i+1)-y) = 0
        % Naj bo H = 1 + a + c in Delta = y - x
        H = 1 + a + c;
        
        % To vodi do kvadratne enačbe A*Delta^2 + B*Delta + C_const = 0
        A = H;
        B = -(H * (delta_i + delta_i1) + b + d_val);
        C_const = fx * delta_i * delta_i1;
        
        % Stabilno reševanje kvadratne enačbe
        diskriminanta = sqrt(B^2 - 4*A*C_const);
        if B > 0
            Delta = -2*C_const / (B + diskriminanta);
        else
            Delta = -2*C_const / (B - diskriminanta);
        end
        
        % Posodobitev približka
        x = x + Delta;
    end
end

function x = racionalni_newton_robni(d, u, ro, k, tip, x0)
    % za vsak slučaj še min pa max lastni vrednosti (AI je zastonj skuhu,
    % tko da zakaj pa ne)
    % Če je ro > 0 deluje samo max, če je ro < 0 deluje samo min
    % v teh primerih je lastna vrednost na prvem ali zadnjem intervalu in 
    % uporabi navadnega racionalnega newtona

    % d: urejen padajoče (d1 > d2 > ... > dn)
    % tip: 'max' za največjo (nad d1) ali 'min' za najmanjšo (pod dn)
    
    n = length(d);
    
    % 1. Določitev začetnega približka in meja, če niso podani
    if strcmp(tip, 'max')
        if nargin < 6 || isempty(x0)
            % Zgornja meja po Gershgorinu ali preko sledi
            x = d(1) + ro * (u' * u); 
        else
            x = x0;
        end
        robni_pol = d(1);
    else % 'min'
        if nargin < 6 || isempty(x0)
            x = d(n) + ro * (u' * u);
        else
            x = x0;
        end
        robni_pol = d(n);
    end

    for step = 1:k
        % Izračun f(x) in f'(x)
        razlike = d - x;
        prispevki = (ro * u.^2) ./ razlike;
        prispevki_der = (ro * u.^2) ./ (razlike.^2);
        
        fx = 1 + sum(prispevki);
        dfx = sum(prispevki_der);
        
        % Za robne intervale vzamemo le en pol (tisti najbližji)
        % f(y) approx (1 + a) + b / (robni_pol - y)
        delta = robni_pol - x;
        b = dfx * delta^2;
        a = (fx - 1) - b/delta;
        
        % Rešujemo: 1 + a + b/(robni_pol - y) = 0
        % y = robni_pol + b / (1 + a)
        x = robni_pol + b / (1 + a);
    end
end

x1 = racionalni_newton(d, u, p, 1, 1);
disp(x1);

% 3 koraki se še vidi razliko, od 4 naprej ostane isto => damo 20 da bo zih
x_majkemi = racionalni_newton(d, u, p, 20, 1);
disp(x_majkemi);

z = u ./ (d - x1);


v_orig = zeros(size(u));
v_orig(idx) = z; 

v_final = v_orig / norm(v_orig, 2);
disp(abs(v_final(1)))
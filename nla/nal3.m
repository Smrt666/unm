A = [1.2094    0.4998    0.9372    0.9284    0.7304    0.4561    1.2465    1.4024    1.1405    0.7302; 
0.4998    1.5025    1.0958    0.8576    0.7749    1.0531    1.4098    0.6539    0.5700    1.3023; 
0.9372    1.0958    0.5086    1.6451    0.2975    1.7235    0.4338    1.1585    1.1617    1.2381; 
0.9284    0.8576    1.6451    1.1705    1.0805    1.2284    0.5148    0.7618    0.9356    0.7314; 
0.7304    0.7749    0.2975    1.0805    1.5583    1.4625    1.0432    1.3437    0.7332    0.5252; 
0.4561    1.0531    1.7235    1.2284    1.4625    0.3313    0.7544    1.0803    0.7996    1.0910; 
1.2465    1.4098    0.4338    0.5148    1.0432    0.7544    1.9016    1.4070    1.1322    0.1541; 
1.4024    0.6539    1.1585    0.7618    1.3437    1.0803    1.4070    0.1689    1.2691    0.4998; 
1.1405    0.5700    1.1617    0.9356    0.7332    0.7996    1.1322    1.2691    1.1594    0.6732; 
0.7302    1.3023    1.2381    0.7314    0.5252    1.0910    0.1541    0.4998    0.6732    0.3678];

function A = jacobi_rotacija(A, i, j)
    % Uniči element (i, j)
    % Preveri, če je element že enak nič
    if abs(A(i,j)) < 1e-15
        return;
    end
    
    % 1. Izračun parametrov rotacije (tau, t, c, s)
    tau = (A(j,j) - A(i,i)) / (2 * A(i,j));
    
    if tau >= 0
        t = 1 / (tau + sqrt(1 + tau^2));
    else
        t = -1 / (-tau + sqrt(1 + tau^2));
    end
    
    c = 1 / sqrt(1 + t^2);
    s = c * t;
    
    % 2. Posodobitev diagonalnih elementov
    % Uporabimo stabilnejše formule za A(i,i) in A(j,j)
    Aii_star = A(i,i) - t * A(i,j);
    Ajj_star = A(j,j) + t * A(i,j);
    
    % 3. Posodobitev ostalih vrstic in stolpcev (le i in j)
    for k = 1:size(A, 1)
        if k ~= i && k ~= j
            A_ki = A(k, i);
            A_kj = A(k, j);
            
            % Rotacija stolpcev
            A(k, i) = c * A_ki - s * A_kj;
            A(k, j) = s * A_ki + c * A_kj;
            
            % Zaradi simetrije posodobimo še vrstice
            A(i, k) = A(k, i);
            A(j, k) = A(k, j);
        end
    end
    
    % 4. Zapišemo nove vrednosti na diagonalo in uničimo tarčni element
    A(i,i) = Aii_star;
    A(j,j) = Ajj_star;
    A(i,j) = 0;
    A(j,i) = 0;
end

A2 = jacobi_rotacija(A, 1, 2);
disp(A2(1,1));

function val = off(A)
    % Izračuna normo nediagonalnih elementov
    A_off = A - diag(diag(A));
    val = sqrt(sum(A_off(:).^2));
end

function [row, col, max_val] = najvecji_element(A)
    % Poišče indekse največjega nediagonalnega elementa po absolutni vrednosti
    n = size(A, 1);
    A_temp = abs(A - diag(diag(A))); % Ignoriramo diagonalo
    [max_val, idx] = max(A_temp(:));
    [row, col] = ind2sub([n, n], idx);
end

function A = klasicni_jacobi(A, toleranca)
    % Izvaja rotacije na največjih elementih, dokler vsi niso pod toleranco
    while true
        [i, j, val] = najvecji_element(A);
        
        % Preverimo pogoj: ali so VSI nediagonalni elementi pod 10^-3
        if val < toleranca
            break;
        end
        
        % Izvedemo rotacijo (uporabimo prej definirano funkcijo)
        A = jacobi_rotacija(A, i, j);
    end
end

Ak = klasicni_jacobi(A, 1e-3);
disp(off(Ak));

function A = ciklicni_jacobi(A, toleranca)
    % Izvaja rotacije po vrsti na trikotniku
    while max(max(abs(A - diag(diag(A))))) >= toleranca
        n = size(A, 1);
        for i = 1:n
            for j = i+1:n
                A = jacobi_rotacija(A, i, j);
            end
        end
    end
end

function val = off_abs(A)
    % Izračuna 1 normo nediagonalnih elementov
    A_off = A - diag(diag(A));
    val = sum(sum(abs(A_off)));
end

Ac = ciklicni_jacobi(A, 1e-2);
disp(off_abs(Ac));
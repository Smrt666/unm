function [output] = mat123(n)
    output = diag(1:n) + triu(ones(n), 1) * 2 + diag(-1 * ones(n-1, 1),-1);
end

% disp(mat123(5));

function [A, A2, a2] = matsq(n)
    A = randi(2*n + 1, n) - (n + 1) * ones(n);
    A2 = A ^ 2;
    a2 = A .* A;
end

% [a, b, c] = matsq(5);
% disp(a);
% disp(b);
% disp(c);

function [A] = vecdiv(a, b)
    b(b == 0) = 1;
    A = a' * b .^ -1;
end

% disp(vecdiv([0, 1, 2, -3], [-3, 0, 2, 1]));

function [y] = fprasaj()
    x = input("Stevilo: ");
    if x < 6
        y = 2;
    elseif 6 <= x && x <= 20
        y = x - 4;
    else
        y = -x;
    end
end

% fprasaj()

function [vsota, delne] = vec_sum(v)
    vsota = sum(v);
    delne = zeros(size(v));
    delne(1) = v(1);
    for i = 2:length(v)
        delne(i) = v(i) + delne(i-1);
    end
end

% [v, d] = vec_sum(randi(100, 20, 1));
% disp(v);
% disp(d);
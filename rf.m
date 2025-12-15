f = @(x) 5*cos(x - exp(x)) - x;

a = 0;
b = 3;
tol = 1e-6;
maxit = 100;

if f(a)*f(b) > 0
    error('Funkcija na robovih intervala nima nasprotnih predznakov.');
end

[x, xs, k] = regula_falsi(f, a, b, tol, maxit);

fprintf('Približek ničle: x_%d = %.10f\n', k, x);

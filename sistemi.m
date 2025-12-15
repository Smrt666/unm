F = @(x) [x(1)^5 + x(2)^5 - 5*x(1)*x(2) - 1; x(1)^2 + x(2)*exp(x(2)^2) - 1];
J = @(v) [
    5*v(1)^4 - 5*v(2),      5*v(2)^4 - 5*v(1);
    2*v(1),                exp(v(2)^2)*(1 + 2*v(2)^2) ];
G = @(x) [nthroot(1 + 5*x(1)*x(2) - x(2)^5, 5); (1 - x(1)^2) / exp(x(2)^2)];

x0 = [0.5; -0.5];
%[x, a, k] = jac_iter(G, x0, 1e-5, 100);
%disp(x);
%disp(k);

%[x, xs, k] = newton(F, J, x0, 1e-5, 100);
%disp(x);
%disp(xs);
%disp(k);

[x, xs, k] = broyden(F, J(x0), x0, 1e-5, 100);
disp(x);
disp(xs);
disp(k);
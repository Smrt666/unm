gs = {@(x) x*x - 2, @(x) sqrt(x + 2), @(x) 1 + 2/x, @(x) (x*x + 2) / (2*x - 1)};

tol = 1e-8;
n = 100;
x0s = linspace(-2, 4, 601);

for j = (1:length(gs))
    g = gs{j};
    ks = zeros(1, 600);
    i = 1;
    for x0 = x0s
        [xr, xarr, k] = iteracija(g, x0, tol, n);
        ks(i) = k;
        i = i + 1;
    end
    figure
    plot(x0s, ks);
    ylim([0, 100]);
end

function [p, dp] = poli0(roots, a)
   p = poly(roots) * a;
   dp = polyder(p);
   x = linspace(min(roots) - 1, max(roots) + 1);
   plot(x, polyval(p, x));
   hold on;
   plot(x, polyval(dp, x));
   hold off;
end

[p, dp] = poli0([-1,2,3], -0.5);
disp(p);
disp(dp);
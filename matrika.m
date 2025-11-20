function [A] = matrika(n)
    A = diag(-(2:2:(2*n))) + diag(n - (1:(n-1)), 1) + diag(n - (1:(n-1)), -1);
end
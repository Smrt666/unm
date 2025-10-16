function [mat,prod1] = Matrika1(n, m, varargin)
    narginchk(1, 2);
    if nargin == 1
        m = n;
    end
    mat = ((1:n)' * (1:m)) ./ (repmat(1:m, n, 1) + repmat(1:n, m, 1)');
    prod1 = prod(mat(:, 1));
end
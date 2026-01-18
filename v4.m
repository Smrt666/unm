function [A] = amat(b)
    A = eye(10) * b;
    for i = 1:10
        for j = 1:10
            A(i, j) = A(i, j) + mod(j - i, 10);
        end
    end
end

A = amat(17);
[q, r] = givens(A);
disp(max(max(abs(q))));
disp(max(max(abs(r))));

% qmax = zeros(0, 101);
% rmax = zeros(0, 101);
% for b = -50:50
%     A = amat(b);
%     [Q, R] = givens(A);
%     qmax(b + 51) = max(max(abs(Q)));
%     rmax(b + 51) = max(max(abs(R)));
% end
% 
% disp(qmax);
% disp(rmax);
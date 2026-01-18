A = [4 1 -1 0;
    1 3 -1 0;
    -1 -1 5 2;
    0 0 2 4];
zac_pribl = ones(4,1);
max_kor = 10;
tol = 1e-16;

[kappa] = obcutljivost(A, zac_pribl, max_kor, tol)

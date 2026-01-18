function[kappa] = obcutljivost(A, zac_pribl, max_kor, tol)
B = A'*A;
[lastna_vr_1, ~, ~] = potencna(B, zac_pribl, tol, max_kor);
[lastna_vr_n, ~, ~] = potencnainv(B, zac_pribl, tol, max_kor);
kappa = sqrt(lastna_vr_1/lastna_vr_n)
end
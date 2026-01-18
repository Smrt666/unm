function [lastna_vr, lastni_vek, st_kor] = potencnainv(A, zac_pribl, tol, max_kor)
B = inv(A);
lastni_vek = zac_pribl;
st_kor = 0;
lastna_vr = lastni_vek'*B*lastni_vek/(lastni_vek'*lastni_vek);
while (norm(B*lastni_vek - lastna_vr*lastni_vek) >= tol) & (st_kor < max_kor)
    st_kor = st_kor + 1;
    y = B*lastni_vek;
    lastni_vek = y/norm(y);
    lastna_vr = lastni_vek'*B*lastni_vek/(lastni_vek'*lastni_vek);
end
lastna_vr = 1/lastna_vr
end
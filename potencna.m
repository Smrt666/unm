function [lastna_vr, lastni_vek, st_kor] = potencna(A, zac_pribl, tol, max_kor)
lastni_vek = zac_pribl;
st_kor = 0;
lastna_vr = lastni_vek'*A*lastni_vek/(lastni_vek'*lastni_vek);
while (norm(A*lastni_vek - lastna_vr*lastni_vek) >= tol) & (st_kor < max_kor)
    st_kor = st_kor + 1;
    y = A*lastni_vek;
    lastni_vek = y/norm(y);
    lastna_vr = lastni_vek'*A*lastni_vek/(lastni_vek'*lastni_vek);
end
end



function A=rowhouse(A,v)

% ROWHOUSE(A,v) za dano m*n matriko A in nenicelen vektor v dimenzije m
% z v(1)=1 mnozi matriko A s P=I-2*v*v'/(v'*v) z leve.

% Bor Plestenjak
% NM 2000/2001

b=-2/(v'*v);
w=b*(A'*v);
A=A+v*w';

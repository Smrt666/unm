function A=rowrot(A,c,s)

% ROWROT(A,c,s) za dano matriko A z dvema vrsticama pomnozi
% matriko A z Givensovo rotacijo [c,s; -s,c]' z leve.

% Bor Plestenjak
% NM 2000/2001

r1=A(1,:); 
r2=A(2,:);
A(1,:)=c*r1-s*r2; 
A(2,:)=s*r1+c*r2;


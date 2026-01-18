function [Q,R]=qrgivensdemo(A)

% QRGIVENSDEMO(A) po korakih naredi QR razcep matrike A
% preko Givensovih rotacij. 
% predpostavljeno je, da je matrika A m x n, kjer je m>=n.

% Bor Plestenjak
% NM 2000/2001

[m,n]=size(A);
Q=eye(m);

for i=1:n
  for j=i+1:m
    [c,s]=givens(A(i,i),A(j,i));
    A([i j],:)=rowrot(A([i j],:),c,s);  
    Q([i j],:)=rowrot(Q([i j],:),c,s);  
    clc
    A
    Q'
    disp('Pritisni enter za nadaljevanje')
    pause
  end
end 

Q=Q';
R=A;

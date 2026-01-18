function [Q,R]=qrhouse(A)

% QRHOUSE(A) naredi QR razcep matrike A preko
% Householderjevih zrcaljenj. V matriki Q vrne vektorje,
% ki predstavljajo Householderjeva zrcaljenja.
% Predpostavljeno je, da je matrika A reda m x n,
% kjer je m>=n.

% Bor Plestenjak
% NM 2000/2001

[m,n]=size(A);
V=zeros(m);

if m>n 
  konec=n;
else
  konec=n-1;
end

for k=1:konec
  v=house(A(k:m,k));
  A(k:m,k:n)=rowhouse(A(k:m,k:n),v);
  V(k+1:m,k)=v(2:m-k+1); % shranimo vse razen 1 v spodnji del A ( sedaj V)
end

R=A;

Q=eye(m);
for j=konec:-1:1
  v(j)=1;
  v(j+1:m)=V(j+1:m,j);
  Q(j:m,j:m)=rowhouse(Q(j:m,j:m),v(j:m));
end
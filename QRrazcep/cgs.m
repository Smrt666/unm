function [Q,R]=cgs(A)

% [Q,R]=CGS(A) vrne QR razcep matrike A preko
% klasicne Gram-Schmidtove ortogonalizacije

% Bor Plestenjak
% NM 2000/2001

[m,n]=size(A);

for i=1:n
   Q(:,i)=A(:,i);
   for j=1:i-1
      R(j,i)=Q(:,j)'*A(:,i);
      Q(:,i)=Q(:,i)-R(j,i)*Q(:,j);
   end
   R(i,i)=norm(Q(:,i));
   Q(:,i)=Q(:,i)/R(i,i);
end

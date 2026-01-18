function [c,s]=givens(a,b)

% GIVENS(a,b) vrne c in s tako da matrika [c s;-s c]' unici drugo komponento
% vektorja [a b]'

% Bor Plestenjak
% NM 2000/2001

if b==0
   c=1; s=0;
else
   if abs(b)>abs(a)
      r=-a/b; s=1/sqrt(1+r*r); c=s*r;
   else
      r=-b/a; c=1/sqrt(1+r*r); s=c*r;
   end
end

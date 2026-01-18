function v=house(x)

% HOUSE(x) za n-dimenzionalni vektor x izracuna n-dimenzionalni vektor v,
% tako da I-2*v*v'/(v'*v) unici vse razen prve komponente v vektorju x.

% Bor Plestenjak
% NM 2000/2001

n=length(x);
u=norm(x);
v=x;
if u~=0
  b=x(1)+signum(x(1))*u;
  v(2:n)=v(2:n)/b;
end
  v(1)=1;

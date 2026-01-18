function s=signum(x)

% SIGNUM(x) je signum funkcija, ki jo potrebujemo pri Householderjevih
% zrcaljenjih. To pomeni SIGNUM(x)=1 za x>=0, sicer pa -1.

% Bor Plestenjak
% NM 2000/2001

if x>=0
  s=1;
else
  s=-1;
end

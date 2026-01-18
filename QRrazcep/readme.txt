Programi za raèunanje QR razcepa matrike A:

[Q,R]=cgs(A) : klasièna Gram-Schmidtova ortogonalizacija

[Q,R]=mgs(A) : modificirana Gram-Schmidtova ortogonalizacija

[Q,R]=QRGivens(A) : QR razcep preko Givensovih rotacij, pomozne funkcije:

   [c,s]=givens(a,b) : izraèuna c in s, ki unièi drugo komponento v vektorju [a;b]
   A=rowrot(A,c,s)   : naredi Givensovo rotacijo z c in s na matriki A, ki ima dve vrstici

[Q,R]=QRHouse(A) : QR razcep preko Householderjevih zrcaljenj, pomozne funkcije:

   v=house(x)        : izracuna v, ki bo prezcralil x v k*e1 
   s=signum(x)       : predznak x (1 za x>=0)
   A=rowhouse(A,v)   : matriko A pomnozi z leve z zrcaljenjem doloèenim z v 

[Q,R]=QRGivensDemo(A) : QR razcep preko Givensovih rotacij po posameznih korakih

[Q,R]=QRHouseDemo(A)  : QR razcep preko Householderjevih zrcaljenj po korakih


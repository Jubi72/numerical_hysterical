% Blatt 01

%Beispielmatrizen aus Beispiel 2.3, A=LR
L= [1, 0, 0;
    2/3, 1, 0;
    1/3, 2, 1];
R= [3, 1, 6;
    0, 1/3, -1;
    0, 0, 1];
b= [2; 7; 4];

%Loesen Sie hier die Gleichung Ax=b mit Hilfe der beiden Funktionen
%vorwaerts und rueckwaerts


x=?


%Implementieren Sie hier die Funktion vorwaerts(L,b), die eine Gleichung
%der Form Lx=b loest. L normierte untere linke Dreiecksmatrix
function [x]=vorwaerts(L,b)
n=size(L,2);
x=zeros(n,1);




end

%Implementieren Sie hier die Funktion rueckwaerts(L,b), die eine Gleichung
%der Form Rx=b loest. R rechte obere Dreiecksmatrix
function [x]=rueckwaerts(R,b)



end
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

y = vorwaerts(L,b)
x = rueckwaerts(R,y)

% Wir haben gesehen, dass es auch eine direkte Summenformel ohne for-Schleife
% gibt, allerdings erst nach erfolgter Implementierung. Wir wollten jetzt dafuer
% nicht den bereits funktionierenden Algorithmus nachtraeglich veraendern.

%Implementieren Sie hier die Funktion vorwaerts(L,b), die eine Gleichung
%der Form Lx=b loest. L normierte untere linke Dreiecksmatrix
function [x]=vorwaerts(L,b)
n=size(L,2);
x=zeros(n,1);
for j = 1:n
    x(j) =  b(j);
    for k = 1:j-1
        x(j) = x(j) - L(j,k) * x(k);
    end
end
end

%Implementieren Sie hier die Funktion rueckwaerts(L,b), die eine Gleichung
%der Form Rx=b loest. R rechte obere Dreiecksmatrix
function [x]=rueckwaerts(R,b)
n = size (R,2);
x = zeros(n,1);
x(n) = b(n)/R(n,n);
for j = n-1:-1:1
    x(j) = b(j);
    for k = j+1:n
        x(j) = x(j) - R(j,k) * x(k);
    end
    x(j) = x(j)/R(j,j);
end
end

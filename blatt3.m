% Blatt 03

% Hilbertmatrizen koennen mit der Funktion hilb erstellt werden
n=8;      H=hilb(n);


%Hinweise: 
% 1) y=noisyVector(b,p)   kann zum Generieren gestoerter rechter Seiten
% verwendet werden
% 2) die Funktionen vorwaerts und rueckwaerts vom Blatt 1 koennen ebenfalls
% verwendet werden
% 3) semilogy  koennen Sie zum plotten verwenden
% 4) Hilfreich koennen auch cond und norm sein
% Beachten Sie bei allen Matlabfunktionen die optionalen Parameter


Nmax=10; 
kond=zeros(Nmax,1);
fehlerH=zeros(Nmax,1);
fehlerX=zeros(Nmax,1); fehlerB=fehlerH; 






%Implementieren Sie hier den Cholesky-Algorithmus
function [L]=cholesky(A)



end





% Funktion zum Addidieren von p% Rauschen 
function y=noisyVector(b,p)
% y=b+noise, wobei noise normierter Fehlervektor (normalverteilt) ist
% der relative Fehler ist ||y-b||/||b|| = p/100

    delta=randn(size(b));                            % Zufallsvektor delta
    noise=delta/norm(delta,'inf')*norm(b,'inf');     % normieren des Fehlervektors   
    y = b+ noise* p/100;                             % additiver Fehler (p%)
end



% Implementation der Vorwaerts- und Rueckwaertselimination 
function [x]=rueckwaerts(R,b)
n=size(R,2);
x=zeros(n,1);
for i=n:-1:1
    x(i)=(b(i)-R(i,i+1:n)*x(i+1:n))/R(i,i);
end
end

% Im Vergleich zu Blatt 1 ist die vorwaerts Funktion leicht modifiziert
% worden fuer allgemeine untere Dreiecksmatrizen L 
function [x]=vorwaerts(L,b)
n=size(L,2);
x=zeros(n,1);
for i=1:n
    x(i)=(b(i)- L(i,1:i-1)*x(1:i-1))/L(i,i);
end
end


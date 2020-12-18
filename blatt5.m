%Blatt 05
eps=1e-8; % Toleranz
kmax=1e8; % maximale Iterationsanzahl (können Sie auch innerhalb der Fkt definieren)
n=51;     % Anzahl der Gitterpunkte
 
% Vektor der gesuchten Unbekannten
u=zeros(n,1);
% Matrix A und rechte Seite f
A=[];
f=[];
% Relaxationsparameter
omega=[];
    

% Hinwise: 
% - Versuchen Sie, mit möglichst wenig Speicher auszukommen (also mit
% sparsen Matrizen zu arbeiten) und eine moeglichst schnellen Algorithmus
% zu implementieren
% - Sie koennen selbstverstaendlich die Unterfunktionen veraendern oder
% zusaetzliche hinzufuegen wenn nötig






% Funktion für die Gauss Seidel Iteration 
% Wenn Sie mehr Eingangs- oder Ausgangsparameter benoetigen, konnen Sie die
% Funktion entsprechend anpassen
function [x]= gaussseidel(A,x,b)



end

 
% Gauss Seidel mit Relaxation
function [x]= SOR(A,x,b,omega)



end


% vorwaerts Funktion vom Blatt 3 fuer allgemeine untere Dreiecksmatrizen L 
% darf selbstverstaendlich ebenfalls angepasst werden
function [x]=vorwaerts(L,b)
n=size(L,2);
x=zeros(n,1);

for i=1:n
    x(i)=(b(i)- L(i,1:i-1)*x(1:i-1))/L(i,i);
end
end



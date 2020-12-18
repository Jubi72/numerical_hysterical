%Blatt 6

global eps;
global maxk;
eps=1e-8;
maxk=50;   % fuer Teilaufgabe 2 

% Loesung
x_sol=[sqrt(2),2];

% Menge an Spalten-Startvektoren fuer das Newton-Verfahren ( Teilaufgabe 2)
X=[20,5,1,-1;20,3,1,-3]; 

% Startpunkt fuer Teilaufgabe 3 % 4
X=[20,20];    
    




%%  Ab hier ist Platz fuer Ihre Funktionen



function x=Newton(x0,func    ) 
% Die Funkion  [x, erg]=Newton(x0,fun, varargin) berechnet Nullstellen der
global maxk;
global eps;
x = x0;
Deltax = 1;
%lastx = 1; % brauchen wir gar nicht
for k = 1:maxk
    if norm(func(x) < eps)
        break
    end
    % ableitung = abl berechnen
    % Delta xk berechnen mit func(x) + abl(x)Deltax = 0
    x = x + Deltax;
end
end 


function []=fkt(x) 
% Funktion aus der Aufgabe 6.3 
% x1^2 + x2^2 = 6 und x2=x1^2


end

%% LR- und Elminationsfunktionen von den verherigen Uebungsblaettern
function [A] = LR(A)
n=size(A,1);
for k=1:n-1
    if abs(A(k,k))<=10^(-6)
        disp("Fehler!")
        break
    else
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
    end
end
end
function [x]=rueckwaerts(R,b)
n=size(R,2);
x=zeros(n,1);
for i=n:-1:1
    x(i)=(b(i)-R(i,i+1:n)*x(i+1:n))/R(i,i);
end
end
function [x]=vorwaerts(L,b)
n=size(L,2);
x=zeros(n,1);
for i=1:n
    x(i)=(b(i)- L(i,1:i-1)*x(1:i-1))/L(i,i);
end
end
%Blatt 05
global eps;
eps=1e-8; % Toleranz
global kmax;
kmax=1e8; % maximale Iterationsanzahl (können Sie auch innerhalb der Fkt definieren)
global n;
n=100;     % Anzahl der Gitterpunkte
global ITER
ITER = zeros(20,1);

h=1/n;    % Gitterweite

% Vektor der gesuchten Unbekannten
u=zeros(n,1);

% Matrix A und rechte Seite f


e = ones(n,1);
A = 1/(h*h*h*h)*spdiags([e -4*e 6*e -4*e e], -2:2,n,n);

for i=1:2               % u_0=u_1=u_(n−1)=u_n= 0
    for j=1:4
        A(i,j)=0;
        A(j,i)=0;
        A(n-i+1,n-j+1)=0;
        A(n-j+1,n-i+1)=0;
    end
end


f=zeros(n,1);           % f = -1 für x in [0.5,0.7]; sonst 0
for i=1:n
    if i/n >= 0.5
        if i/n <= 0.7
            f(i) = -1;
        end
    end
end

g = gaussseidel(A,u,f)

% Relaxationsparameter
omega=1:19;

for i=1:size(omega,2)
    s = SOR(A,u,f,omega(i));
end

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
L=tril(A,-1);
R=triu(A,1);
D=diag(diag(A));
global eps;
global kmax;
global n;
tic;
for k=1:kmax
    lastx = x;
    y = -inv(L+D)*R*x+(L+D)\b;
    x(3:n-2) = y(3:n-2);
    if norm(lastx-x) <= eps*norm(x)
        k
        break
    end
end
toc
end

 
% Gauss Seidel mit Relaxation
function [x]= SOR(A,x,b,omega)
global ITER;
kmax = 1000;
global eps;
L=tril(A,-1);
R=triu(A,1);
D=diag(diag(A));
omegaindex = omega;
omega = omega/10; %omega wird nur als index uebergeben, hier berechnen
global n;         %wir den eigentlichen relaxationsparameter
for k=1:kmax
    lastx = x;
    y = -inv (omega*L+D) * ( (1-omega) *R) * x + (omega*L+D) \ b;
    x(3:n-2) = y(3:n-2);
    if norm(lastx-x) <= eps*norm(x)
        ITER(omegaindex)=k;
        break
    end
end
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



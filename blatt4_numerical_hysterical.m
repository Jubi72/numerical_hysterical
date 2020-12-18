%Blatt 04

kmax=6;                         % bestimmt die maximale Groesse von A 
global Nmax;
Nmax = 20000;                 % maximale Anzahl von Iterationen fuer Ihre Schleife
global eps;
eps=1e-8;
global ITER;
ITER = zeros(kmax,3);
global OFFS;
OFFS = zeros(kmax,3);
global k;
for k= 1:kmax
    m=2^k;  n=m^2;
    % Matrix A als Blockmatrix (vom Blatt 2)
    e = ones(m,1);
    B = spdiags([-e 4*e -e], -1:1, m, m);   % Bandmatrix B
    E=speye(m);                             % Identitaet                     
    D=spdiags([-e -e], [-1,1],m,m);
    A=kron(E,B) + kron(D,E);                % Matrix A als Kroneckerprodukt
   
    xb=ones(n,1);
    
    b=A*xb;                                 % Rechte Seite des LGS
   
    % Hinweise
    % statt einer for-schleife  kann auch eine while-schleife verwendet werden
    % Setzten Sie sicherheitshalber eine maximale Anzahl Nmax maximaler
    % Iterationen fest
    % Zum Schaetzer fuer q -  setzen Sie ein sinnvolles q_0 fest
    
    x = zeros(n,1);
    format long;
    Jacobi(A,x,b);
end
tiledlayout(2,1) 
nexttile
semilogy(ITER)
legend('nach Residuum','Differenz','Fehlerschätzung','Location','northwest')
title("Iterationen")
nexttile
semilogy(OFFS)
legend('nach Residuum','Differenz','Fehlerschätzung','Location','northwest')
title("Fehler")
saveas(1, 'ergebnis_blatt4', 'png')
ITER
OFFS
%% Implementation des Jacobi-Verfahrens

function [x] = Jacobi(A,x,b)
global Nmax;
global eps;
global k;
global ITER;
global OFFS;
fertig = zeros(3,1);
r0 = b - A * x;
res = eps*norm(r0);
lastx = x;
for n = 0:Nmax
    beforelastx = lastx;
    lastx = x;
    r = b - A * x;
    deltax = linsolve(full(diag(A)),r);
    x = x + deltax;
    q = norm(x - lastx)/norm(lastx - beforelastx);
    if norm(r) <= res                              % nach Residuum (a)
        if fertig(1) == 0
            ITER(k,1) = n;
            OFFS(k,1) = 1 - x(1);
            fertig(1) = 1;
        end
    end
    if norm(x-lastx) <= eps*norm(x)                % nach Differenz zweier ben. Iter.
        if fertig(2) == 0
            ITER(k,2) = n;
            OFFS(k,2) = 1 - x(1);
            fertig(2) = 1;
        end
    end
    if norm(x - lastx) <= eps*(1-q)*norm(lastx)    % nach Fehlerschätzung
        if fertig(3) == 0
            ITER(k,3) = n;
            OFFS(k,3) = 1 - x(1);
            fertig(3) = 1;
        end
    end
    if fertig == ones(3,1)
        break
    end
end
end
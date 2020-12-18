%Blatt 04

kmax=6;                     % bestimmt die maximale Groesse von A 
Nmax=20000;                 % maximale Anzahl von Iterationen fuer Ihre Schleife
eps=1e-8;


    k=2;
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
    
    
   
%% Implementation des Jacobi-Verfahrens

function [x] = Jacobi(A,x,b)

end
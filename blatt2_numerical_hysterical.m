
kmax=5;  % Hier sollten Sie bei einem ersten Test zunächst mit einem k< 7 starten
% Lasttests wurden auf anderem Rechner durchgefuehrt

% Vektoren zeit und fehler für die Teilaufgabe 3
% Vektor fehler umbenannt in n, da dort ns statt Fehler eingetragen werden
zeit=zeros(kmax,1); ns=zeros(kmax,1);  

for k = 1:kmax
    m=2^k;  n=m^2;
    
    % Blockmatrizen bauen mit Hilfe von Kroneckerprodukten:
    % sparse Matrizen  (Tipp: mit full(B) können Sie wieder
    % die volle Matrix erhalten um Sie zu prüfen/visualisieren)
    e = ones(m,1);
    B = spdiags([-e 4*e -e], -1:1, m, m);   % Bandmatrix B
    E=speye(m);                             % Identität                     
    D=spdiags([-e -e], [-1,1],m,m);
    A=kron(E,B) + kron(D,E);                % Matrix A als Kroneckerprodukt
    
    % Nun können Sie starten mit Ihrer LR-Zerlegung
    % Hilfreich können folgende Befehler sein:
    % tic  und toc für die Zeitmessung 
    % tril und triu  für Dreiecksmatrizen
    
    tic;
    LR(A);
    zeit(k) = toc; 
    % Wir messen wirklich nur die Zeit fuer die LR-Zerlegung, nicht mehr
    ns(k) = n;
end

%% Graphische Ausgabe: Hier können Sie Ihr Ergebnis zeichen und speichern
% Speichen  des Fenters 1 können Sie zum Beispiel mit Hilfe von
% saveas(1, 'ergebnis_blatt2', 'png')

plot(ns, zeit);
saveas(1, 'ergebnis_blatt2', 'png');
% s. o.
    
%% Funktion LR zur Berechnung der LR-Zerleung der Matrix A
function [A] = LR(A)
n = size(A,1);
L = eye(n);
R = A;
for j = 1:n-1
    if A(j,j) == 0
        fprintf('Fehler: Diagonale der Eingabe enthaelt 0\n')
        return
    end
    R(j,j+1:n) = A(j,j+1:n);
    % folgende Zeile taucht zwar im Skript im Algorithmus 2 auf,
    % ist aber offensichtlich algorithmisch ueberfluessig
    % A(j,j+1:n) = R(j,j+1:n);
    L(j+1:n,j) = A(j+1:n,j) / A(j,j);
    A(j+1:n,j) = L(j+1:n,j);
    A(j+1:n,j+1:n) = A(j+1:n,j+1:n) - L(j+1:n,j) * R(j,j+1:n);
end
end

function [f, g] = myobjective(z, A)
% myobjective.m: Calcola f_A(z) e g_A(z) per la matrice A e un punto complesso z.
%
% Input:
%   z - Punto (numero complesso) per il quale calcolare f_A(z) e g_A(z).
%   A - Matrice quadrata (reale o complessa).
%
% Output:
%   f - Valore della funzione f_A(z) = det(A - zI).
%   g - Valore di g_A(z) definito come 1/trace((A - zI)^{-1}).
%
% Descrizione:
%   La funzione calcola il determinante della matrice B = A - zI utilizzando
%   la fattorizzazione LU con pivoting parziale e il calcolo del determinante
%   della matrice di permutazione, secondo quanto indicato nella consegna:
%       f_A(z) = det(P) * prod(diag(U))
%   Inoltre, calcola l'inversa di B sfruttando la fattorizzazione LU per ottenere
%   g_A(z) = 1/trace((A - zI)^{-1}).
%

% Determina la dimensione n della matrice A
n = size(A,1);

% Costruisce la matrice B = A - zI.
% Questo sistema (A - zI)x = 0 è alla base della definizione di autovalore.
B = A - z * eye(n);

% Esegue la fattorizzazione LU con pivoting parziale in forma vettoriale.
% [L, U, piv] = lu(B, 'vector') restituisce:
%   L: matrice triangolare inferiore,
%   U: matrice triangolare superiore,
%   piv: vettore che codifica la matrice di permutazione P.
[L, U, piv] = lu(B, 'vector');

% Calcola il determinante del fattore di permutazione P.
% Il vettore 'piv' rappresenta una permutazione delle righe.
% Per ottenere det(P) si conta il numero di scambi (swaps) necessari per 
% riportare il vettore alla sequenza ordinata, in quanto:
%   det(P) = (-1)^(numero di scambi)
num_swaps = 0;
temp = piv; 
for i = 1:n
    while temp(i) ~= i
        j = temp(i);
        % Scambia gli elementi per "riordinare" il vettore
        [temp(i), temp(j)] = deal(temp(j), temp(i));
        num_swaps = num_swaps + 1;
    end
end
detP = (-1)^num_swaps;

% Calcola il determinante di U come il prodotto degli elementi sulla sua diagonale.
% Poiché U è triangolare, det(U) = prod(diag(U)).
detU = prod(diag(U));

% Calcola f_A(z) = det(B) = det(P)*det(U).
% Questo segue dalla proprietà della fattorizzazione LU: P*B = L*U.
f = detP * detU;

% Per calcolare g_A(z), occorre calcolare l'inversa di B.
% Ricostruiamo la matrice di permutazione P a partire dal vettore 'piv'.
Pmat = eye(n);      % Inizializza la matrice identità
Pmat = Pmat(piv,:);  % Riordina le righe secondo il vettore 'piv' per ottenere P

% Calcola l'inversa di B sfruttando i fattori della LU.
% Dato che B = Pmat'*L*U, si ha:
%   B^{-1} = U^{-1} * L^{-1} * Pmat
% Per evitare l'uso esplicito di inv(), si risolve il sistema lineare.
Binv = U \ (L \ Pmat);

% Calcola la traccia di B^{-1}, cioè la somma degli elementi diagonali.
theTrace = trace(Binv);

% Calcola g_A(z) come g = 1/trace(B^{-1}).
% Questo segue direttamente dalla definizione data nella consegna:
%   g_A(z) = 1/trace((A - zI)^{-1})
g = 1 / theTrace;

end

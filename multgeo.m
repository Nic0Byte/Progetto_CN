function k = multgeo(A, l, toll)
% multgeo.m: Calcola la molteplicità geometrica di un autovalore l di A.
%
% Input:
%   A    - Matrice quadrata (reale o complessa).
%   l    - Autovalore (complesso) per cui si vuole calcolare la molteplicità geometrica.
%   toll - Tolleranza usata per decidere se un elemento diagonale di U è considerato nullo.
%
% Output:
%   k    - Molteplicità geometrica, ovvero il numero di righe "nulle" ottenute
%          dalla fattorizzazione LU di B = A - lI.
%
% Riferimenti dalla consegna:
%   - La molteplicità geometrica di un autovalore l è definita come:
%         mult_geo(l) = n - rank(A - lI)
%     In altre parole, è il numero di vettori linearmente indipendenti che
%     risolvono (A - lI)x = 0.
%
%   - Si suggerisce di utilizzare la fattorizzazione LU (con pivoting parziale)
%     per determinare la presenza di righe "nulle". In particolare, gli elementi
%     sulla diagonale della matrice U, se inferiori ad una certa soglia (toll),
%     sono considerati numericamente nulli, indicando una riduzione del rango.
%
% Procedura:
%   1. Estrae la dimensione n della matrice quadrata A.
%   2. Costruisce la matrice B = A - lI.
%   3. Calcola la fattorizzazione LU di B. La matrice U contiene i pivot sulla diagonale.
%   4. Estrae il vettore dei pivot (valori sulla diagonale di U) e ne calcola il modulo.
%   5. Conta quanti pivot hanno modulo minore di toll. Questo conteggio è la molteplicità
%      geometrica dell'autovalore l.

    % Estrae il numero di righe (dimensione) della matrice A.
    n = size(A, 1);

    % Costruisce la matrice B = A - lI.
    % Questo corrisponde al sistema (A - lI)x = 0, per cui una soluzione non nulla
    % indica che l è un autovalore di A.
    B = A - l * eye(n);
    
    % Esegue la fattorizzazione LU con pivoting parziale della matrice B.
    % La funzione lu restituisce L, U e la matrice di permutazione, ma qui ci interessa
    % solo U, in cui i pivot sono posizionati lungo la diagonale.
    [~, U, ~] = lu(B);
    
    % Estrae il vettore dei valori sulla diagonale di U e ne calcola il valore assoluto.
    % Secondo la consegna, i pivot con valore assoluto inferiore alla soglia toll
    % sono considerati come indicatori di righe "nulle" nella matrice, e quindi
    % contribuiscono alla molteplicità geometrica dell'autovalore.
    diagU = abs(diag(U));
    
    % Conta quanti elementi della diagonale di U hanno modulo inferiore a toll.
    % Questo numero equivale alla dimensione dello spazio nullo di (A - lI), ovvero:
    %   mult_geo(l) = n - rank(A - lI) ~ numero di pivot piccoli.
    k = sum(diagU < toll);
end

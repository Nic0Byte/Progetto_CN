function [lambda, m, flag] = multalg(A, l0, toll, it, maxit)
% multalg.m: Calcola un autovalore (lambda) e la sua molteplicità algebrica (m)
%            per la matrice A usando il metodo di Newton standard e modificato.
%
% Input:
%   A     - Matrice quadrata reale o complessa.
%   l0    - Punto di partenza per il metodo di Newton.
%   toll  - Tolleranza per il test d'arresto su |z_{k+1} - z_k|.
%   it    - Numero di iterazioni per il metodo di Newton standard.
%   maxit - Numero massimo di iterazioni per il metodo di Newton modificato.
%
% Output:
%   lambda - Autovalore approssimato.
%   m      - Molteplicità algebrica stimata.
%   flag   - Flag di successo: 1 se la convergenza è stata raggiunta, 0 altrimenti.
%
% Descrizione:
%   1. Si applica il metodo di Newton standard per 'it' iterazioni.
%      Se la differenza tra due iterazioni è minore di 'toll', si assume 
%      convergenza con molteplicità 1.
%
%   2. Se non si raggiunge la convergenza, si eseguono due ulteriori iterazioni 
%      per calcolare gli step s1 e s2 e stimare il rapporto degli step,
%      che approssima (m-1)/m. Da questo si ottiene un valore iniziale per la molteplicità.
%
%   3. Viene quindi applicato il metodo di Newton modificato, dove lo step 
%      viene moltiplicato per una ipotesi di molteplicità (mm). Se la convergenza 
%      viene raggiunta entro 'maxit' iterazioni (o se il numero totale di chiamate 
%      a myobjective supera 10*maxit), si restituiscono i valori trovati; altrimenti, 
%      si incrementa mm e si ripete il processo.
%
% Nota: La funzione myobjective(z, A) restituisce f_A(z) e g_A(z), dove:
%       g_A(z) = 1/trace((A - zI)^{-1}), come specificato nella consegna.
%
% Riferimenti dalla consegna:
%   - La relazione asintotica: |s_k| / |s_{k+1}| ~ (m-1)/m.
%   - Se dopo 10*maxit chiamate a myobjective non si raggiunge la convergenza,
%     si imposta flag = 0.

% Inizializzazione del punto di partenza per il metodo di Newton standard
z = l0;
[~, g_cur] = myobjective(z, A);  % Utilizziamo solo g_A(z), ignorando f_A(z)

% ==============================
% Metodo di Newton Standard
% ==============================
for k = 1:it
    % Calcola il nuovo punto con l'aggiornamento del metodo di Newton
    z_new = z + g_cur;
    
    % Se lo step è inferiore alla tolleranza, consideriamo la convergenza
    if norm(z_new - z) < toll
        lambda = z_new; % Convergenza raggiunta
        m = 1;          % Radice semplice: molteplicità 1
        flag = 1;
        return;
    end
    
    % Aggiorna il punto e ricalcola g_A(z)
    z = z_new;
    [~, g_cur] = myobjective(z, A);
end

% =====================================================
% Stima della molteplicità tramite rapporto tra step
% =====================================================
% Dopo 'it' iterazioni, se la convergenza non è stata raggiunta,
% eseguiamo due ulteriori iterazioni per ottenere gli step s1 e s2.
z1 = z;
[~, g1] = myobjective(z1, A);
z2 = z1 + g1;

% Se converge in questo passaggio, assumiamo molteplicità 1
if norm(z2 - z1) < toll
    lambda = z2;
    m = 1;
    flag = 1;
    return;
end

[~, g2] = myobjective(z2, A);
z3 = z2 + g2;

% Calcola gli step s1 e s2
s1 = z2 - z1;
s2 = z3 - z2;
ratio = norm(s1) / norm(s2);   % Approssimazione: |s_k|/|s_{k+1}| ~ (m-1)/m

% Controllo per evitare divisione per zero: se ratio è praticamente 1,
% il denominatore della formula m = 1/(1 - ratio) sarebbe quasi zero.
if abs(ratio - 1) < 1e-14
    % Se ratio è molto vicino a 1, assegniamo un valore predefinito a m_guess
    m_guess = 5;
else
    % Calcola m_guess dalla relazione (m-1)/m = ratio, ovvero m = 1/(1 - ratio)
    m_guess = 1 / (1 - ratio);
end

% Imposta m_start come il valore intero di m_guess, garantendo m_start >= 1
m_start = floor(m_guess);
if m_start < 1
    m_start = 1;
end

% =====================================================
% Metodo di Newton Modificato per molteplicità > 1
% =====================================================
% Limita il numero totale di chiamate a myobjective per evitare cicli infiniti
maxCalls = 10 * maxit;
numCalls = 0;

% Ciclo esterno: prova per diverse ipotesi di molteplicità (mm) da m_start fino a 20
for mm = m_start : 20
    % Ripristina il punto di partenza per ogni nuovo tentativo
    zcurr = l0;
    [~, gcurr] = myobjective(zcurr, A);
    numCalls = numCalls + 1;
    
    % Ciclo interno: esegue fino a maxit iterazioni del metodo di Newton modificato
    for iterM = 1:maxit
        % Aggiornamento: z_{k+1} = z_k + mm * g_A(z_k)
        znext = zcurr + mm * gcurr;
        
        % Se lo step è inferiore alla tolleranza, la convergenza è stata raggiunta
        if norm(znext - zcurr) < toll
            lambda = znext;
            m = mm;  % La molteplicità stimata è mm
            flag = 1;
            return;
        end
        
        % Aggiorna il punto corrente
        zcurr = znext;
        [~, gcurr] = myobjective(zcurr, A);
        numCalls = numCalls + 1;
        
        % Se il numero totale di chiamate supera il limite, esci con fallimento
        if numCalls > maxCalls
            lambda = zcurr;
            m = mm;
            flag = 0;
            return;
        end
    end
end

% Se non viene raggiunta la convergenza entro i cicli previsti,
% restituisce gli ultimi valori calcolati e imposta flag = 0.
lambda = zcurr; % Ultimo valore calcolato
m = mm;         % Ultima molteplicità testata
flag = 0;
end

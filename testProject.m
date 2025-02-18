%% TEST 1: Autovalori Ben Spaziati (Molteplicità 1)
% Creiamo J come matrice diagonale (tutti blocchi 1x1) con autovalori ben spaziati.
lambda_vals = [-5, -2, 0, 3, 7];
J1 = [];
for i = 1:length(lambda_vals)
    % Ogni blocco è di dimensione 1x1
    J1 = blkdiag(J1, lambda_vals(i));
end
n1 = size(J1,1);
[Q1, ~] = qr(randn(n1));
A1 = Q1' * J1 * Q1;

% Visualizziamo gli autovalori
eig_A1 = eig(A1);
figure;
plot(real(eig_A1), imag(eig_A1), 'bo', 'MarkerSize',10, 'LineWidth',2);
title('Test 1: Autovalori Ben Spaziati (Molteplicità 1)');
xlabel('Parte Reale'); ylabel('Parte Immaginaria'); grid on;
disp('Test 1: Autovalori Ben Spaziati:');
disp(eig_A1);

% Calcoliamo le molteplicità geometriche attese (dovrebbero essere tutte 1)
for i = 1:length(lambda_vals)
    mg = multgeo(A1, lambda_vals(i), 1e-8);
    disp(['Autovalore ', num2str(lambda_vals(i)), ' -> Molteplicità geometrica: ', num2str(mg)]);
end

%% TEST 2: Autovalori Clusterizzati (Molteplicità Bassa: 1-2)
% Creiamo una matrice con un cluster attorno ad un autovalore, ad es. 7.
% Per simulare il clustering usiamo:
%   - Un blocco di Jordan 2x2 per l'autovalore 7 (molteplicità algebrica 2, geometrica 1)
%   - Altri blocchi 1x1 per valori vicini a 7 (es. 7.05 e 6.95)
J2 = [];
J2 = blkdiag(J2, [7, 1; 0, 7]);   % Blocco di Jordan 2x2 per 7
J2 = blkdiag(J2, 7.05);            % Blocco 1x1 per 7.05
J2 = blkdiag(J2, 6.95);            % Blocco 1x1 per 6.95
n2 = size(J2,1);
[Q2, ~] = qr(randn(n2));
A2 = Q2' * J2 * Q2;

% Visualizziamo gli autovalori
eig_A2 = eig(A2);
figure;
plot(real(eig_A2), imag(eig_A2), 'ro', 'MarkerSize',10, 'LineWidth',2);
title('Test 2: Autovalori Clusterizzati attorno a 7 (Molteplicità 1-2)');
xlabel('Parte Reale'); ylabel('Parte Immaginaria'); grid on;
disp('Test 2: Autovalori Clusterizzati:');
disp(eig_A2);

% Molteplicità geometrica per 7 (dovrebbe essere 1, a causa del blocco 2x2)
mg_cluster = multgeo(A2, 7, 1e-8);
disp(['Autovalore 7 (clusterizzato) -> Molteplicità geometrica: ', num2str(mg_cluster)]);

% Stima della molteplicità algebrica tramite multalg (dovrebbe stimare 2 per il blocco 2x2)
it = 5; maxit = 20;
[lambda_est2, m_est2, flag2] = multalg(A2, 7, 1e-8, it, maxit);
disp(['Test 2: multalg -> lambda stimato: ', num2str(lambda_est2), ...
      ', Molteplicità algebrica stimata: ', num2str(m_est2), ', Flag: ', num2str(flag2)]);

%% TEST 3: Autovalori con Molteplicità Alta (3-7)
% Creiamo una matrice J con blocchi di Jordan di dimensione maggiore:
%   - Un autovalore 3 semplice (blocco 1x1)
%   - Un autovalore 9 con blocco 4x4 (molteplicità algebrica 4, geometrica 1)
%   - Un autovalore 5 con blocco 2x2 (molteplicità algebrica 2, geometrica 1)
J3 = [];
J3 = blkdiag(J3, 3);   % Blocco 1x1 per 3
% Blocco di Jordan 4x4 per 9:
J_block9 = 9 * eye(4) + diag(ones(3,1), 1);
J3 = blkdiag(J3, J_block9);
% Blocco di Jordan 2x2 per 5:
J_block5 = 5 * eye(2) + diag(ones(1,1), 1);
J3 = blkdiag(J3, J_block5);
n3 = size(J3,1);
[Q3, ~] = qr(randn(n3));
A3 = Q3' * J3 * Q3;

% Visualizziamo gli autovalori
eig_A3 = eig(A3);
figure;
plot(real(eig_A3), imag(eig_A3), 'ks', 'MarkerSize',10, 'LineWidth',2);
title('Test 3: Autovalori con Molteplicità Alta (3-7)');
xlabel('Parte Reale'); ylabel('Parte Immaginaria'); grid on;
disp('Test 3: Autovalori con Molteplicità Alta:');
disp(eig_A3);

% Calcoliamo le molteplicità geometriche:
mg_3 = multgeo(A3, 3, 1e-8);
mg_9 = multgeo(A3, 9, 1e-8);
mg_5 = multgeo(A3, 5, 1e-8);
disp(['Autovalore 3 -> Molteplicità geometrica: ', num2str(mg_3)]);
disp(['Autovalore 9 -> Molteplicità geometrica: ', num2str(mg_9)]);
disp(['Autovalore 5 -> Molteplicità geometrica: ', num2str(mg_5)]);

% Stima della molteplicità algebrica per l'autovalore 9 tramite multalg.
% Partiamo da un valore iniziale vicino a 9 (es. 9.1).
[lambda_est3, m_est3, flag3] = multalg(A3, 9.1, 1e-8, 5, 20);
disp(['Test 3: multalg per autovalore 9 -> lambda stimato: ', num2str(lambda_est3), ...
      ', Molteplicità algebrica stimata: ', num2str(m_est3), ', Flag: ', num2str(flag3)]);

%% COMMENTI GENERALI SUL TEST
%
% TEST 1:
% - Viene creata una matrice J1 diagonale con autovalori ben spaziati: -5, -2, 0, 3, 7.
% - La matrice A1 viene ottenuta tramite una trasformazione ortogonale (A1 = Q1' * J1 * Q1),
%   che preserva gli autovalori.
% - Vengono calcolati e visualizzati gli autovalori di A1.
% - Per ciascun autovalore, si calcola la molteplicità geometrica usando la funzione multgeo,
%   che dovrebbe restituire 1 per ciascun autovalore, in quanto ogni blocco è 1x1.
%
% TEST 2:
% - Viene creata una matrice J2 che simula un clustering attorno all'autovalore 7.
% - In J2, si include un blocco di Jordan 2x2 per l'autovalore 7 (molteplicità algebrica 2, 
%   ma molteplicità geometrica 1) e due blocchi 1x1 per valori leggermente diversi (7.05 e 6.95).
% - La matrice A2 viene ottenuta come A2 = Q2' * J2 * Q2.
% - Vengono visualizzati gli autovalori e si calcola la molteplicità geometrica per l'autovalore 7,
%   che dovrebbe essere 1 a causa della presenza del blocco 2x2.
% - La funzione multalg viene utilizzata per stimare la molteplicità algebrica di 7, che dovrebbe essere 2.
%
% TEST 3:
% - Viene creata una matrice J3 con blocchi di Jordan di dimensioni maggiori:
%   - Un blocco 1x1 per l'autovalore 3.
%   - Un blocco 4x4 per l'autovalore 9, che implica una molteplicità algebrica 4 ma geometrica 1.
%   - Un blocco 2x2 per l'autovalore 5, che implica una molteplicità algebrica 2 ma geometrica 1.
% - La matrice A3 è ottenuta tramite una trasformazione ortogonale (A3 = Q3' * J3 * Q3),
%   mantenendo gli autovalori invariati.
% - Gli autovalori di A3 vengono visualizzati, e per ciascuno vengono calcolate le molteplicità geometriche.
% - La funzione multalg viene utilizzata per stimare la molteplicità algebrica dell'autovalore 9,
%   partendo da un valore iniziale vicino a 9 (9.1), e la stima dovrebbe essere 4.
%
% In tutti i test, le trasformazioni ortogonali (QR) garantiscono che la struttura spettrale
% delle matrici J venga preservata in A. Inoltre, le funzioni multgeo e multalg sono verificate
% confrontando i risultati ottenuti con le aspettative teoriche (definite dalle dimensioni dei blocchi
% di Jordan costruiti).

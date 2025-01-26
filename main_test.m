% main_test.m
clear; close all; clc;

% 1) Intervallo e funzione da interpolare
a = -1; 
b =  1;
f = @(x) 1 ./ (x - 1.69);   % Esempio come da testo

% 2) Parametri di test
gradi = [2, 4, 8, 16];      % Esempio
z = linspace(a, b, 1001)';  % Griglia fitta per errori
fz = f(z);

% Preallocazione vettori per risultati
err_equi = zeros(size(gradi));
err_leja  = zeros(size(gradi));
leb_equi  = zeros(size(gradi));
leb_leja  = zeros(size(gradi));

% 3) Loop sui vari gradi
for k = 1:length(gradi)
    n = gradi(k);
    
    % 3.1) Nodi equispaziati
    x_equi = linspace(a, b, n+1)';
    
    % 3.2) Nodi Leja (Algoritmo 1: DLP)
    % prima creo un vettore fitto X su cui estraggo i (n+1) nodi di Leja
    X = linspace(a, b, 500)';   % 500 punti "candidati"
    x_leja = DLP(X, n);         % Algoritmo DLP
    
    % 4) Costruisco la Vandermonde classica per i nodi equi:
    Ve = buildVandermonde(x_equi);
    fe = f(x_equi);
    ce = Ve \ fe;  % risoluzione
    
    % 5) Valuto l’interpolante su z
    p_equi = zeros(size(z));
    for j=1:n+1
        p_equi = p_equi + ce(j)*z.^(j-1);
    end
    
    % 6) Ripeto per i nodi di Leja
    Vl = buildVandermonde(x_leja);
    fl = f(x_leja);
    cl = Vl \ fl;
    p_leja = zeros(size(z));
    for j=1:n+1
        p_leja = p_leja + cl(j)*z.^(j-1);
    end
    
    % 7) Calcolo errori
    err_equi(k) = max(abs(fz - p_equi));
    err_leja(k) = max(abs(fz - p_leja));
    
    % 8) Calcolo costante di Lebesgue
    leb_equi(k) = leb_con(z, x_equi);
    leb_leja(k) = leb_con(z, x_leja);
end

% 9) Visualizzo o stampo i risultati
disp('   n     err_equi    err_leja    leb_equi    leb_leja');
for k=1:length(gradi)
    fprintf('%4d   %1.3e    %1.3e    %1.3e    %1.3e\n', ...
        gradi(k), err_equi(k), err_leja(k), leb_equi(k), leb_leja(k));
end

% 10) Eventuali plot
figure;
semilogy(gradi, err_equi, 'o-', gradi, err_leja, 's-');
legend('Equispaziati','Leja','Location','NorthWest');
xlabel('Grado n'); ylabel('Errore max in scala log');
title('Confronto errori di interpolazione');

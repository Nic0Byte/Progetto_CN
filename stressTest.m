%% STRESS TEST: Matrice con elevata molteplicità e cluster
J_stress = [];
J1 = 1*eye(7) + diag(ones(6,1),1);
J_stress = blkdiag(J_stress, J1);
J2 = 1.001*eye(5) + diag(ones(4,1),1);
J_stress = blkdiag(J_stress, J2);
J_stress = blkdiag(J_stress, 2);
n_stress = size(J_stress,1);
[Q_stress, ~] = qr(randn(n_stress));
A_stress = Q_stress' * J_stress * Q_stress;
eig_stress = eig(A_stress);
figure;
plot(real(eig_stress), imag(eig_stress), 'ms', 'MarkerSize',10, 'LineWidth',2);
title('Autovalori - Stress Test');
xlabel('Parte Reale'); ylabel('Parte Immaginaria');
grid on;

tol = 1e-10;
mg1 = multgeo(A_stress, 1, tol);
mg1_001 = multgeo(A_stress, 1.001, tol);
mg2 = multgeo(A_stress, 2, tol);
fprintf('Molteplicità geometrica per lambda=1 (attesa = 1): %d\n', mg1);
fprintf('Molteplicità geometrica per lambda=1.001 (attesa = 1): %d\n', mg1_001);
fprintf('Molteplicità geometrica per lambda=2 (attesa = 1): %d\n', mg2);

it = 5; maxit = 50;
[lambda_est1, m_est1, flag1] = multalg(A_stress, 1.01, tol, it, maxit);
[lambda_est2, m_est2, flag2] = multalg(A_stress, 1.001, tol, it, maxit);
[lambda_est3, m_est3, flag3] = multalg(A_stress, 2.1, tol, it, maxit);
fprintf('Risultati multalg per lambda vicino a 1: lambda = %g, molteplicità = %d, flag = %d\n', lambda_est1, m_est1, flag1);
fprintf('Risultati multalg per lambda vicino a 1.001: lambda = %g, molteplicità = %d, flag = %d\n', lambda_est2, m_est2, flag2);
fprintf('Risultati multalg per lambda vicino a 2: lambda = %g, molteplicità = %d, flag = %d\n', lambda_est3, m_est3, flag3);

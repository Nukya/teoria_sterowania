%% Zadanie 6 - Sterowanie z minimalnym wydatkiem (L1)
% Dodatkowe ograniczenia:
%   - narzucone wartości u_0 = u_{N-1} = 0,
%   - ograniczone zmiany sterowania: |u_{k+1} - u_k| <= 0.02.

clear; clc; close all;

%% Ustawienia LaTeX
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Parametry modelu
m     = 1.0;
beta  = 0.1;
Delta = 0.5;
tf    = 10.0;
N     = round(tf/Delta);

x0 = [0; 0];
xf = [1; 0];

%% Dyskretyzacja ZOH
[A,B,C] = zoh_discretize(Delta, m, beta);

%% Warunek końcowy na stan
R   = reachability(A,B,N);
rhs = xf - (A^N)*x0;

%% Zmienna decyzyjna LP: z = [u; v]
f = [zeros(N,1); ones(N,1)];

% Podstawowe równanie: R*u = rhs  -> [R, 0]*[u; v] = rhs
Aeq = [R, zeros(2,N)];
beq = rhs;

% Dodatkowe równania: u_0 = 0, u_{N-1} = 0
e1      = zeros(1,N); e1(1)   = 1;
eN      = zeros(1,N); eN(end) = 1;
Aeq_add = [e1, zeros(1,N);
           eN, zeros(1,N)];
beq_add = [0; 0];

Aeq = [Aeq; Aeq_add];
beq = [beq; beq_add];

%% Ograniczenia |u_k| <= v_k oraz v_k >= 0
I = eye(N);
Z = zeros(N);
Aineq_abs = [ I, -I;
             -I, -I;
              Z, -I];
bineq_abs = zeros(3*N,1);

%% Ograniczenia na zmiany sterowania: |u_{k+1} - u_k| <= 0.02
D = zeros(N-1,N);
for k = 1:N-1
    D(k,k)   = -1;
    D(k,k+1) =  1;
end
Aineq_slope = [ D, zeros(N-1,N);
               -D, zeros(N-1,N)];
bineq_slope = 0.02*ones(2*(N-1),1);

Aineq = [Aineq_abs;
         Aineq_slope];
bineq = [bineq_abs;
         bineq_slope];

%% Rozwiązanie LP
opts = optimoptions('linprog','Display','none');
z = linprog(f, Aineq, bineq, Aeq, beq, [], [], opts);
u = z(1:N);

%% Symulacja układu
[xd, yd] = simulate_discrete(A,B,C,u,x0); %#ok<NASGU>
Ac = [0 1;
      0 -beta/m];
Bc = [0;
      1/m];
[t_c, x_c, u_c] = simulate_continuous_zoh(Ac,Bc,u,Delta,x0);

t_u = (0:N-1)*Delta;
t_x = (0:N)*Delta;

%% Wykresy
figure;
subplot(3,1,1);
stairs(t_u, u, 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('$u_k$');
title('Zadanie 6: sterowanie dyskretne $u(k)$');
grid on;

subplot(3,1,2);
plot(t_c, x_c(:,1), 'LineWidth', 1.5); hold on;
plot(t_c, x_c(:,2), 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('$x(t)$');
legend('$x_1(t)$','$x_2(t)$','Location','best');
title('Zadanie 6: przebiegi stanu w czasie ciągłym');
grid on;

subplot(3,1,3);
plot(x_c(:,1), x_c(:,2), 'LineWidth', 1.5);
xlabel('$x_1$'); ylabel('$x_2$');
title('Zadanie 6: trajektoria w przestrzeni stanu');
grid on;

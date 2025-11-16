%% Zadanie 5 - Sterowanie z minimalnym wydatkiem (norma L1), bez dodatkowych ograniczeń
% W tym zadaniu wyznaczamy sterowanie, które minimalizuje sumę modułów u_k
% (tzw. minimal fuel), przy spełnionym warunku końcowym na stanie x_N = x_f.

clear; clc; close all;

%% Ustawienia LaTeX dla wykresów
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Parametry modelu
m     = 1.0;        % masa wózka [kg]
beta  = 0.1;        % współczynnik tarcia lepkiego
Delta = 0.5;        % okres próbkowania [s]
tf    = 10.0;       % czas symulacji [s]
N     = round(tf/Delta);  % liczba kroków sterowania

x0 = [0; 0];        % stan początkowy
xf = [1; 0];        % wymagany stan końcowy po czasie tf

%% Dyskretyzacja ZOH
[A,B,C] = zoh_discretize(Delta, m, beta);

%% Budowa macierzy osiągalności i warunku końcowego
R   = reachability(A,B,N);
rhs = xf - (A^N)*x0;   % prawa strona równania x_N = A^N x0 + R u

%% Formułowanie problemu LP dla normy L1
% Zmienna decyzyjna:
%   z = [u; v], gdzie:
%     u ∈ R^N – sterowanie w kolejnych krokach,
%     v ∈ R^N – zmienne pomocnicze spełniające v_k >= |u_k|.
%
% Funkcja celu: min sum_k v_k  ->  f = [0; ...; 0; 1; ...; 1].
f = [zeros(N,1); ones(N,1)];

% Ograniczenie równościowe: R*u = rhs
% W przestrzeni z = [u; v] zapisujemy jako [R, 0] * [u; v] = rhs.
Aeq = [R, zeros(2,N)];
beq = rhs;

% Ograniczenia nierównościowe zapewniające v_k >= |u_k| oraz v_k >= 0:
%   u_k - v_k <= 0    -> [ I  -I ]
%  -u_k - v_k <= 0    -> [-I  -I ]
%        -v_k <= 0    -> [ 0  -I ]
I = eye(N);
Z = zeros(N);
Aineq = [ I, -I;
         -I, -I;
          Z, -I];
bineq = zeros(3*N,1);

%% Rozwiązanie zadania LP za pomocą linprog
opts = optimoptions('linprog','Display','none');
z = linprog(f, Aineq, bineq, Aeq, beq, [], [], opts);
u = z(1:N);   % wyodrębniamy sterowanie u z wektora z

%% Symulacja układu dyskretnego
[xd, yd] = simulate_discrete(A,B,C,u,x0); 

%% Symulacja układu ciągłego (dla potrzeb rysunków)
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
title('Zadanie 5: sterowanie dyskretne $u(k)$');
grid on;

subplot(3,1,2);
plot(t_c, x_c(:,1), 'LineWidth', 1.5); hold on;
plot(t_c, x_c(:,2), 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('$x(t)$');
legend('$x_1(t)$','$x_2(t)$','Location','best');
title('Zadanie 5: przebiegi stanu w czasie ciągłym');
grid on;

subplot(3,1,3);
plot(x_c(:,1), x_c(:,2), 'LineWidth', 1.5);
xlabel('$x_1$'); ylabel('$x_2$');
title('Zadanie 5: trajektoria w przestrzeni stanu');
grid on;

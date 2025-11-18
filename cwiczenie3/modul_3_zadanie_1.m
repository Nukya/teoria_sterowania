close all
clear
clc

%% --- Globalne ustawienia LaTeX + FontSize ---
nfonts = 14;        % czcionka dla osi, legend
nfontslatex = 18;   % czcionka dla opisów, tytułów

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

set(groot,'defaultAxesFontSize',nfonts);
set(groot,'defaultTextFontSize',nfontslatex);
set(groot,'defaultLegendFontSize',nfonts);

%% --- Wczytanie danych ---
load plantData01.mat      % zawiera F, G, H, Q, R
load signalsData01.mat    % zawiera U, Y, X_true

%% Macierze modelu
A = F;
B = G;
C = H;

%% --- Rozpakowanie struktur z sygnałami ---
U = U.u(:);          % 120×1
Y = Y.y(:);          % 121×1
xTrue = X_true.x;    % 3×121

%% --- Długości sygnałów ---
N_u = length(U);     % 120
N_y = length(Y);     % 121
N = N_u;             % filtr działa dla 120 próbek

%% --- Parametry układu ---
n = size(A,1);

%% --- Inicjalizacja filtru Kalmana ---
x_hat_plus = zeros(n,1);
P_plus = 10 * eye(n);

Xminus = zeros(n,N);
Xplus  = zeros(n,N);

%% --- Filtr Kalmana ---
for k = 1:N

    % Predykcja
    x_hat_minus = A * x_hat_plus + B * U(k);
    P_minus = A * P_plus * A' + Q;

    % Korekcja
    K = P_minus * C' / (C * P_minus * C' + R);
    x_hat_plus = x_hat_minus + K * (Y(k) - C * x_hat_minus);
    P_plus = (eye(n) - K*C) * P_minus * (eye(n) - K*C)' + K*R*K';

    % Zapis
    Xminus(:,k) = x_hat_minus;
    Xplus(:,k)  = x_hat_plus;
end

%% --- Wykresy ---
figure('units','normalized','outerposition',[0 0 0.6 1])

%% Panel (a) – sygnały u(k), y(k)
subplot(4,1,1)
plot(0:N_y-1, Y, 'r'); hold on
plot(0:N_u-1, U, 'b');
grid on
xlabel('$k$')
ylabel('$u(k),\; y(k)$')
legend('$y(k)$','$u(k)$')

%% Panel (b) – x1
subplot(4,1,2)
plot(0:N_y-1, xTrue(1,:), 'k'); hold on
plot(0:N-1, Xplus(1,:), 'b');
plot(0:N-1, Xminus(1,:), 'r');
grid on
xlabel('$k$')
ylabel('$x_{1}$')
legend('$x_1$','$\hat{x}_1^+$','$\hat{x}_1^-$')

%% Panel (c) – x2
subplot(4,1,3)
plot(0:N_y-1, xTrue(2,:), 'k'); hold on
plot(0:N-1, Xplus(2,:), 'b');
plot(0:N-1, Xminus(2,:), 'r');
grid on
xlabel('$k$')
ylabel('$x_{2}$')
legend('$x_2$','$\hat{x}_2^+$','$\hat{x}_2^-$')

%% Panel (d) – x3
subplot(4,1,4)
plot(0:N_y-1, xTrue(3,:), 'k'); hold on
plot(0:N-1, Xplus(3,:), 'b');
plot(0:N-1, Xminus(3,:), 'r');
grid on
xlabel('$k$')
ylabel('$x_{3}$')
legend('$x_3$','$\hat{x}_3^+$','$\hat{x}_3^-$')

clear
close all
clc

set(groot, 'defaultTextInterpreter', 'latex'); 
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');

%% Parametry obiektu
m = 1;
beta = 0.1;
Delta = 0.5;

alpha = -beta/m;
b = 1/m;

%% Macierze układu dyskretnego A,B
A = [1, Delta;
     0, exp(alpha*Delta)];

B = [ (Delta^2)/2;
      (exp(alpha*Delta)-1)/alpha ];

%% Dane zadania
x0 = [0; 0];
xf = [1; 0];

tf = 10;
N = tf / Delta;

%% Budowa macierzy R wg wzoru (5)
R = zeros(2, N);
for j = 1:N
    R(:, j) = A^(N-j) * B;
end

%% ZADANIE 7 — minimalny wydatek sterowania (norma L1)
% przy ograniczeniach:
%   u0 = u_{N-1} = 0
%   u_{k+1} - u_k <= 0.02

cvx_begin quiet
    variable mu(N)     % sterowanie u0 ... u_{N-1}
    variable t(N)      % zmienne pomocnicze do linearyzacji normy L1

    minimize( sum(t) ) % minimalizacja sumy wartości bezwzględnych sterowania

    subject to
        % warunek osiągnięcia stanu końcowego
        R * mu == xf - A^N * x0;

        % u0 = 0 oraz u_{N-1} = 0
        mu(1) == 0;
        mu(N) == 0;

        % ograniczenie różnicowe: u_{k+1} - u_k <= 0.02
        for k = 1:N-1
            mu(k+1) - mu(k) <= 0.02;
        end

        % L1 linearizacja
        -t <= mu <= t;
        t >= 0;
cvx_end


%% Symulacja układu dyskretnego (równanie stanu)
Xd = zeros(2, N+1);
Xd(:,1) = x0;

for k = 1:N
    Xd(:,k+1) = A*Xd(:,k) + B*mu(k);
end

%% Symulacja układu ciągłego (ode45)
Ac = [0 1; 0 alpha];
Bc = [0; b];

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[tPlot, Xc] = ode45(@(t,x) Ac*x + Bc*mu(min(floor(t/Delta)+1, N)), ...
                    [0 tf], x0, opts);

u_cont = arrayfun(@(t) mu(min(floor(t/Delta)+1, N)), tPlot);

%% --- Rysunki ---
% Sterowanie
figure;
subplot(2,1,1)
t_disc = 0:Delta:tf-Delta;
plot(t_disc, mu, 'k.', 'MarkerSize', 18);
grid on
xlabel('$t\,[s]$')
ylabel('$u(t)$')
ylim([-0.2 0.4])
title('Sygnał sterowania u(t) - układ z czasem dyskretnym')


subplot(2,1,2)
plot(tPlot, u_cont, 'LineWidth', 2)
grid on
xlabel('$t\,[s]$')
ylabel('$u(t)$')
ylim([-0.2 0.4])
title('Sygnał sterowania u(t) - układ z czasem ciągłym')

% zmienne stanu
figure;
subplot(2,1,1)
t_disc = 0:Delta:tf;

plot(t_disc, Xd(1,:), '.', 'MarkerSize', 18, 'Color', [0 0 1]);  % x1 — niebieskie kropki
hold on
plot(t_disc, Xd(2,:), '.', 'MarkerSize', 18, 'Color', [1 0 0]);  % x2 — czerwone kropki

grid on
xlabel('$t\,[s]$')
ylabel('$x(t)$')
ylim([0 1])
legend('$x_1$','$x_2$')
title('Zmienne stanu - układ z czasem dyskretnym')


subplot(2,1,2)
plot(tPlot, Xc(:,1), 'LineWidth', 2)
hold on
plot(tPlot, Xc(:,2), 'LineWidth', 2)
grid on
xlabel('$t\,[s]$')
ylabel('$x(t)$')
ylim([0 1])
legend('$x_1$','$x_2$')
title('Zmienne stanu - układ z czasem ciągłym')

% Trajektoria w przestrzeni stanu
figure;
subplot(1,2,1)
plot(Xd(1,:), Xd(2,:), '.', 'MarkerSize', 18);
grid on
xlabel('$x_1$')
ylabel('$x_2$')
ylim([0 1])
xlim([0 1])
title('Trajektoria w przestrzeni stanu - układ z czasem dyskretnym')

subplot(1,2,2)
plot(Xc(:,1), Xc(:,2), 'LineWidth',2)
grid on
xlabel('$x_1$')
ylabel('$x_2$')
ylim([0 1])
xlim([0 1])
title('Trajektoria w przestrzeni stanu - układ z czasem ciągłym')

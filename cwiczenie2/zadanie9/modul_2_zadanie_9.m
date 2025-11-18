clear
close all
clc

set(groot, 'defaultTextInterpreter', 'latex'); 
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');

%% Parametry obiektu
m = 1;
beta = 0.1;
Delta = 0.1;

alpha = -beta/m;
b = 1/m;

%% Macierze układu dyskretnego A,B
A = [1, Delta;
     0, exp(alpha*Delta)];

B = [ (Delta^2)/2;
      (exp(alpha*Delta)-1)/alpha ];

C = [1 0];   % wyjście y = x1

%% Dane zadania
x0 = [0; 0];

kf = 100;
tf = kf * Delta;

omega = 2*pi/10;

%% Trajektoria referencyjna (84)
% y_ref(k) = sin(k * omega * Delta), k = 1,...,kf
yref = sin((1:kf)' * omega * Delta);

%% Budowa macierzy S i R (wzór (25): y = S x0 + R mu)

S = zeros(kf, 2);
for k = 1:kf
    S(k,:) = C * A^k;          % wiersz odpowiada y_k
end

R = zeros(kf, kf);
for i = 1:kf
    for j = 1:i
        R(i,j) = C * A^(i-j) * B;  % wpływ u_{j-1} na y_i
    end
end

t_disc = (1:kf) * Delta;   % chwile próbkowania dla y i e
t_u    = (0:kf-1) * Delta; % chwile próbkowania dla u

%% ============================
%  ZADANIE 9 — część 1: gamma = 0
%% ============================

gamma = 0;

mu = (R' * R + gamma * eye(kf)) \ (R' * (yref - S * x0));

%% Symulacja układu dyskretnego
Xd = zeros(2, kf + 1);
Xd(:,1) = x0;
Yd = zeros(kf, 1);

for k = 1:kf
    Xd(:,k+1) = A * Xd(:,k) + B * mu(k);
    Yd(k) = C * Xd(:,k+1);   % y_k = C x_k, tutaj k = 1..kf
end

e = Yd - yref;

%% --- Rysunek 26: gamma = 0 ---
figure;

% (a) sterowanie u(t)
subplot(2,2,1)
plot(t_u, mu, 'k.', 'MarkerSize', 12)
grid on
xlabel('$t\,[s]$')
ylabel('$u(t)$')
title('(a) Sterowanie $u(t)$')

% (b) uchyb e(t)
subplot(2,2,2)
plot(t_disc, e, 'k.', 'MarkerSize', 10)
grid on
xlabel('$t\,[s]$')
ylabel('$e(t)$')
title('(b) Uchyby $e(t) = y(t) - y_{ref}(t)$')

% (c) y_ref(t)
subplot(2,2,3)
plot(t_disc, yref, 'k.', 'MarkerSize', 10)
grid on
xlabel('$t\,[s]$')
ylabel('$y_{ref}(t)$')
title('(c) Odpowiedz referencyjna $y_{ref}(t)$')

% (d) y(t)
subplot(2,2,4)
plot(t_disc, Yd, 'k.', 'MarkerSize', 10)
grid on
xlabel('$t\,[s]$')
ylabel('$y(t)$')
ylim([-1 1])
title('(d) Odpowiedz ukladu $y(t)$ na sterowanie $u(t)$')


%% ============================
%  ZADANIE 9 — część 2: gamma = 0.1
%% ============================

gamma = 0.1;

mu2 = (R' * R + gamma * eye(kf)) \ (R' * (yref - S * x0));

%% Symulacja układu dyskretnego
Xd2 = zeros(2, kf + 1);
Xd2(:,1) = x0;
Yd2 = zeros(kf, 1);

for k = 1:kf
    Xd2(:,k+1) = A * Xd2(:,k) + B * mu2(k);
    Yd2(k) = C * Xd2(:,k+1);
end

e2 = Yd2 - yref;

%% --- Rysunek 27: gamma = 0.1 ---
figure;

% (a) sterowanie u(t)
subplot(2,2,1)
plot(t_u, mu2, 'k.', 'MarkerSize', 12)
grid on
xlabel('$t\,[s]$')
ylabel('$u(t)$')
ylim([-1 4])
title('(a) Sterowanie $u(t)$')

% (b) uchyb e(t)
subplot(2,2,2)
plot(t_disc, e2, 'k.', 'MarkerSize', 10)
grid on
xlabel('$t\,[s]$')
ylabel('$e(t)$')
title('(b) Uchyby $e(t) = y(t) - y_{ref}(t)$')

% (c) y_ref(t)
subplot(2,2,3)
plot(t_disc, yref, 'k.', 'MarkerSize', 10)
grid on
xlabel('$t\,[s]$')
ylabel('$y_{ref}(t)$')
ylim([-1 1])
title('(c) Odpowiedz referencyjna $y_{ref}(t)$')

% (d) y(t)
subplot(2,2,4)
plot(t_disc, Yd2, 'k.', 'MarkerSize', 10)
grid on
xlabel('$t\,[s]$')
ylabel('$y(t)$')
ylim([-1 1])
title('(d) Odpowiedz ukladu $y(t)$ na sterowanie $u(t)$')

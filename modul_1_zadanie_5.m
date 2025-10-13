clear;
clc;
close all;

% Parametry układu 
L = 180e-6;     % [H]
C = 5.4e-6;     % [F]
R = 9;          % [Ohm]
Vin = 15;       % [V]
fs = 40e3;      % [Hz]
ws = 2*pi*fs;   % [rad/s]

% Czas symulacji 
tspan = [0 3e-3]; % [s]

% Funkcja przełączająca 
u = @(t) 0.5*(1 - sign(sin(ws*t)));

% Uproszczone macierze stanu zakładając R1=R2=R3=0
A = [0      -1/L;
    1/C -1/(R*C)];

B = [1/L;
    0];

b = [0; 0];

f = [0; Vin/(C*R)];

% Równanie stanu 
dxdt = @(t, x) (A*x + B*u(t) + f);

% Warunki początkowe 
x0 = [0; 0]; % [i_l, v_C]

% Symulacja 
[t, x] = ode45(dxdt, tspan, x0);


iL = x(:, 1);
vC = x(:, 2); 

% Sekcja kreślenia wykresu 
figure;
subplot(2,1,1);
plot(t, iL, 'LineWidth', 1.2);
xlabel('t [s]');
ylabel('i_L [A]');
title('Przebieg prądu dławika');
grid on;

subplot(2,1,2);
plot(t, vC, 'LineWidth', 1.2);
xlabel('t [s]');
ylabel('v_C [V]');
title('Przebieg napięcia na kondensatorze');
grid on;




clear;
clc;
close all;

% Surowe dane

% L = 180uH
% C = 5.4 uF
% R = 9 Ohm
% Vi = 15 V 
% fs = 40 kHz

% Dane w podstawowych jednostkach SI
L = 180e-6;     % [H]           180uH
C = 5.4e-6;     % [F]           5.4 uF
R = 9;          % [Ohm]         9 Ohm
Vi = 15;        % [V]           15 V 
fs = 40e3;      % [Hz]          40 kHz
ws = 2*pi*fs;   % [rad/s]

% Rozpatrywany czas
tspan = [0 3e-3]; % [s]

% Funkcja przełączania
u = @(t) (1-sign(sin(ws*t)))/2;

% Uproszczone macierze stanu dla R1=R2=R3=0
A = [0      -1/L;
    1/C -1/(R*C)];

B = [1/L;
    0];

b = [0; 0];

f = [0; Vi/(C*R)];

% Równanie stanu 
dxdt = @(t, x) (A*x + B*u(t) + f);

% Warunki początkowe 
x0 = [0; 0]; % [i_l, v_C]

[t, x] = ode45(dxdt, tspan, x0);

iL = x(:, 1);
vC = x(:, 2); 

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




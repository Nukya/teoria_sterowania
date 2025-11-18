clc; 
clear; 
close all;

% Dane
L = 197e-6;          % H
C = 100e-9;          % F
V0 = 1.95;           % V
Vs = 14;             % V
fs = 40e3;           % Hz
ws = 2*pi*fs;        % rad/s

% Równania
f = @(t, x) [
    (1/L)*(-x(2) - V0*sign(x(1)) + Vs*sign(sin(ws*t)));  % di/dt
    (1/C)*x(1)                                           % dv/dt
];

x0 = [0; 0];   % [i(0); v(0)]

% Przedział
tspan = [0 0.003];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, x] = ode45(f, tspan, x0, opts);

iL = x(:, 1);
vC = x(:, 2);

% wykresy
atrybutyTextu = {'Interpreter', 'latex', 'FontSize', 16};

figure;
subplot(2,1,1);
plot(t, iL, 'LineWidth', 1.5);
grid on;
xlabel('$t\,\mathrm{[ms]}$', atrybutyTextu{:});
ylabel('$i_L\,\mathrm{[A]}$', atrybutyTextu{:});
title('Prad dlawika $i_L(t)$', atrybutyTextu{:});


subplot(2,1,2);
plot(t, vC, 'LineWidth', 1.5);
grid on;
xlabel('$t\,\mathrm{[ms]}$', atrybutyTextu{:});
ylabel('$v_C\,\mathrm{[V]}$', atrybutyTextu{:});
title('Napiecie na kondensatorze $v_c(t)$', atrybutyTextu{:});

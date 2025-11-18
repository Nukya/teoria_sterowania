clc; 
clear; 
close all;

% Dane
L = 197e-6;          % H
C = 100e-9;          % F
Vs = 14;             % V
fs = 40e3;           % Hz
ws = 2*pi*fs;        % rad/s
R = 1.6;             % ohm
C0 = 1e-3;           % F

% Równania
f = @(t, x) [ (1/L)*(-x(2) - x(3)*sign(x(1)) + Vs*sign(sin(ws*t)));   % di/dt
              (1/C)*x(1);                                             % dv/dt
              (1/C0)*(abs(x(1)) - x(3)/R) ];                          % dV0/dt

x0 = [0; 0; 0];   % [i(0); v(0), v0(0)]

% Przedział
tspan = [0 0.003];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, x] = ode45(f, tspan, x0, opts);

iL = x(:, 1);
vC = x(:, 2);
v0C = x(:, 3);

% wykresy
atrybutyTextu = {'Interpreter', 'latex', 'FontSize', 16};

figure;
subplot(3,1,1);
plot(t, iL, 'LineWidth', 1.5);
grid on;
xlabel('$t\,\mathrm{[ms]}$', atrybutyTextu{:});
ylabel('$i_L\,\mathrm{[A]}$', atrybutyTextu{:});
title('Prad dlawika $i_L(t)$', atrybutyTextu{:});


subplot(3,1,2);
plot(t, vC, 'LineWidth', 1.5);
grid on;
xlabel('$t\,\mathrm{[ms]}$', atrybutyTextu{:});
ylabel('$v_C\,\mathrm{[V]}$', atrybutyTextu{:});
title('Napiecie na kondensatorze $v_c(t)$', atrybutyTextu{:});

subplot(3,1,3);
plot(t, v0C, 'LineWidth', 1.4); 
grid on;
xlabel('$t\,\mathrm{[ms]}$', atrybutyTextu{:});
ylabel('$V_0\,\mathrm{[V]}$', atrybutyTextu{:});
title('Napiecie wyjsciowe $V_0(t)$', atrybutyTextu{:});


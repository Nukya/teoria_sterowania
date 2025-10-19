clear; 
close all; 
clc;

% Parametr układu
m = 1.0;

% Definicja równań van der Pola
f = @(t, x) [ x(2);
              m*(1 - x(1)^2)*x(2) - x(1) ];

% Siatka dla pola wektorowego
x1_min = -4;  x1_max = 4;
x2_min = -4;  x2_max = 4;
stepSize = 0.4;

[x1, x2] = meshgrid(x1_min:stepSize:x1_max, x2_min:stepSize:x2_max);

u = x2;
v = m*(1 - x1.^2).*x2 - x1;

% Normalizacja wektorów
n = sqrt(u.^2 + v.^2);
u_n = u ./ n;
v_n = v ./ n;

% Rysowanie portretu fazowego
figure;
quiver(x1, x2, u_n, v_n, 0.5, 'k');
xlabel('x_1');
ylabel('x_2');
title('Portret fazowy dla układu van der Pola');
axis equal; grid on; hold on;

tFinal = 40;
options = odeset('RelTol',1e-9,'AbsTol',1e-12);

% Lista warunków początkowych
x0_list = [
    2.0, 0.0;   % Cykl graniczny 
    0.5, 0.0;
   -2.0, 0.0;
    1.5, 1.0;
];

for i = 1:size(x0_list,1)
    x0 = x0_list(i,:)';
    [t, x] = ode45(f, [0 tFinal], x0, options);

    if all(abs(x0 - [2;0]) < 1e-10)
        plot(x(:,1), x(:,2), 'b', 'LineWidth', 2.0);
    else
        plot(x(:,1), x(:,2), 'r', 'LineWidth', 1.0);
    end
end

legend('Pole wektorowe', 'Trajektorie', 'Location', 'best');
axis tight;
clear; 
clc; 
close all;

% Definicja przedziału i kroku
t = linspace(0,5,1e+3);


% Definicja funkcji
vc = @(t) 2*exp(-t) - 2*(1 + t).*exp(-2*t);
i = @(t) -exp(-t) + (1 + 2*t).*exp(-2*t);

figure;
plot(t, vc(t), 'LineWidth', 1.5); 
hold on;
plot(t, i(t), '--', 'LineWidth', 1.5);
grid on;

atrybutyTextu = {'Interpreter', 'latex', 'FontSize', 16};

xlabel('$t$', atrybutyTextu{:})
ylabel('$v_c \, ,\; i$', atrybutyTextu{:});

title('$v_c(t)$ i $i(t)$ na przedziale $[0,5]$', atrybutyTextu{:});

legend({'$v_c(t)$','$i(t)$'}, atrybutyTextu{:});
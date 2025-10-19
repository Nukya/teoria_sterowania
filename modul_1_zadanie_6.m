clear
close all
clc

nfontslatex = 18;
nfonts = 14;
m = 1.0;   % parametr nieliniowości

% Definicja równania różniczkowego
f = @(t,x) [x(2); m*(1 - x(1)^2)*x(2) - x(1)];

% Funkcje składowe dla pola wektorowego
f1 = @(x1,x2) x2;
f2 = @(x1,x2) m*(1 - x1.^2).*x2 - x1;

% Zakres i siatka
x1_min = -4.0;
x1_max = 4.0;
x2_min = -4.0;
x2_max = 4.0;
stepSize = 0.3;

[X1,X2] = meshgrid(x1_min:stepSize:x1_max, x2_min:stepSize:x2_max);
Y1 = arrayfun(f1,X1,X2);
Y2 = arrayfun(f2,X1,X2);

% Normalizacja wektorów (na potrzeby quivera)
L = sqrt(Y1.^2 + Y2.^2);
Y1_n = Y1 ./ L;
Y2_n = Y2 ./ L;

% Rysunek pola wektorowego
figure
quiver(X1,X2,Y1_n,Y2_n,'b')
daspect([1 1 1])
hold on
xlabel('$x_1$','Interpreter','latex','FontSize',nfontslatex)
ylabel('$x_2$','Interpreter','latex','FontSize',nfontslatex)
title('Trajektorie równania Van der Pola','Interpreter','latex','FontSize',nfontslatex)

tInit = 0.0;
tFinal = 20.0;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% kilka różnych warunków początkowych
xInits = [  2   0;
           -2   0;
            0   2;
            0  -2;
            1   1;
           -1  -1];

colors = lines(size(xInits,1));

for i = 1:size(xInits,1)
    xInit = xInits(i,:)';
    [t,X] = ode45(f,[tInit,tFinal],xInit,options);
    plot(X(:,1),X(:,2),'Color',colors(i,:),'LineWidth',2.0)
end

axis tight
legend(arrayfun(@(i) sprintf('x_0 = [%g, %g]', xInits(i,1), xInits(i,2)), 1:size(xInits,1), 'UniformOutput', false), ...
       'Location','bestoutside')
grid on
hold off

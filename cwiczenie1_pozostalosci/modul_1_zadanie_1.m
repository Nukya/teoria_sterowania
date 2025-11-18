clear
close all
clc
nfontslatex = 18;
nfonts = 14;
m = 1.0;
f = @(t,x) [x(2); m*(1-x(1)^2)*x(2)-x(1)];
tInit = 0.0;
tFinal = 20.0;
xInit = [2; 0];
options = odeset('RelTol',1e-14,'AbsTol',1e-14);
[t,X] = ode45(f,[tInit,tFinal],xInit,options);
figure
plot(t,X(:,1),'r','LineWidth',2.0)
hold on
grid on
plot(t,X(:,2),'b','LineWidth',2.0)
set(gca,'FontSize',nfonts)
xlabel('$t$','Interpreter','Latex', ...
'FontSize',nfontslatex)
ylabel('$x_{1},\,x_{2}$', ...
'Interpreter','Latex','FontSize',nfontslatex)
legend({'$x_{1}$','$x_{2}$'}, ....
'Interpreter','Latex', ...
'FontSize',nfontslatex,'Location','Best')
print('VDP.eps','-depsc','-r600')
print('VDP.jpg','-djpeg','-r600')
print('VDP.pdf','-dpdf','-r600')
%% Parametry
m = 1;                    % kg
beta = 0.1;               % Ns/m
Delta = 0.5;              % s (okres próbkowania)
N = 20;                   % horyzont (tf = 10 s)
x0 = [0;0];
xf = [1;0];

%% Model ciagly i dyskretyzacja ZOH (dokładna)
alpha = -beta/m;
Ac = [0 1; 0 alpha];
Bc = [0; 1/m];

% Dokładne A,B dla ZOH przez macierz powiększoną (unikamy ręcznych wzorów)
M  = [Ac Bc; zeros(1,3)];
Md = expm(M*Delta);
A  = Md(1:2,1:2);
B  = Md(1:2,3);
C  = [1 0];

%% Minimalna energia: równanie osiągalności x_N = A^N x0 + R*u
% Macierz osiągalności na horyzoncie N
n = size(A,1);
R = zeros(n,N);
for i = 1:N
    R(:,i) = A^(N-i) * B;
end

d = xf - A^N * x0;

% Rozwiązanie o najmniejszej normie (energia minimalna): u* = R'*(R*R')^{-1}*d
u_min = R' * ((R*R') \ d);

% (opcjonalnie) to samo przez QP: min 1/2 u'Hu s.t. R u = d
H = eye(N); f = zeros(N,1);
opts = optimoptions('quadprog','Display','off');
u_qp = quadprog(H,f,[],[],R,d,[],[],[], opts);

%% Symulacja dyskretna (rekurencyjnie i przez lsim)
x = zeros(n,N+1); x(:,1) = x0;
for k = 1:N
    x(:,k+1) = A*x(:,k) + B*u_min(k);
end

t = (0:N)*Delta;

% (alternatywnie) lsim dla układu dyskretnego – dostaniemy też ślad stanów
sysd = ss(A,B,C,0,Delta);
[~, t_lsim, x_lsim] = lsim(sysd, u_min, (0:N-1)'*Delta, x0); 

%% Wykresy
figure; 
stairs((0:N-1)*Delta, u_min, 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('u_k'); title('Sterowanie minimalnoenergetyczne');

figure;
plot(t, x(1,:), 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('x_1 = p [m]'); title('Położenie wózka');

figure;
plot(t, x(2,:), 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('x_2 = \dot p [m/s]'); title('Prędkość wózka');

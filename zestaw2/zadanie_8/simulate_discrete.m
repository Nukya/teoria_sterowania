function [x,y] = simulate_discrete(A,B,C,u,x0)
    N = numel(u);
    n = size(A,1);
    x = zeros(n, N+1);
    x(:,1) = x0(:);
    for k = 1:N
        x(:,k+1) = A*x(:,k) + B*u(k);
    end
    y = (C*x).';
end

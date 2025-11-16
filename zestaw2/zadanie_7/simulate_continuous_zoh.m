function [t_all, x_all, u_all] = simulate_continuous_zoh(Ac,Bc,mu,Delta,x0)
    N = numel(mu);
    t_all = [];
    x_all = [];
    u_all = [];
    x_curr = x0(:);
    t_curr = 0;
    for k = 1:N
        u_k = mu(k);
        rhs = @(t,x) Ac*x + Bc*u_k;
        [t_seg, x_seg] = ode45(rhs, [t_curr, t_curr + Delta], x_curr);
        t_all = [t_all; t_seg(1:end-1)];
        x_all = [x_all; x_seg(1:end-1,:)];
        u_all = [u_all; u_k * ones(length(t_seg)-1,1)];
        x_curr = x_seg(end,:).';
        t_curr = t_curr + Delta;
    end
    t_all = [t_all; t_curr];
    x_all = [x_all; x_curr.'];
    u_all = [u_all; mu(end)];
end

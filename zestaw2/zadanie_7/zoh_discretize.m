function [A,B,C] = zoh_discretize(Delta, m, beta)
    alpha = -beta/m;
    b = 1/m;
    if abs(alpha) < 1e-12
        A = [1, Delta; 0, 1];
        B = [0.5*b*Delta^2; b*Delta];
    else
        exp_aD = exp(alpha * Delta);
        psi    = (exp_aD - 1) / alpha;
        A = [1, psi; 0, exp_aD];
        B = [((exp_aD - 1)/alpha^2 - Delta/alpha)*b;
              ((exp_aD - 1)/alpha)*b];
    end
    C = [1 0];
end

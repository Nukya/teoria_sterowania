function R = reachability(A,B,N)
    n = size(A,1);
    R = zeros(n,N);
    for j = 1:N
        R(:,j) = (A^(N-j)) * B;
    end
end

function [R, Rp] = gen_r(E, N)
    R = zeros(E*(N+1),E*N+1);
    q = eye(N + 1);
    
    ri = 1:N+1:E*(N+1);
    ci = 1:N:E*N+1;
    ci = ci(1:end-1);
    
    for i = 1:length(ri)
        R(ri(i):ri(i)+N,ci(i):ci(i)+N) = q;
    end
    
    Rp = R(:,1:end-1);
    Rp(end,1) = 1;
end
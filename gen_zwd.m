function [z, w, p_gll, D] = gen_zwd(N)

    % LIBRARY OF COEFFICIENTS

    if N > 9 || N < 3
        msg = 'Enter polynomial range between 3 and 9';
        error(msg)
    end
    
    legrendre_coeff = zeros(10,7);
    legrendre_diff = legrendre_coeff;

    legrendre_coeff(:,1) = 0.5*[0 0 0 0 0 0 5 0 -3 0]; % N = 3
    legrendre_coeff(:,2) = 0.125*[0 0 0 0 0 35 0 -30 0 3];
    legrendre_coeff(:,3) = 0.125*[0 0 0 0 63 0 -70 0 15 0];
    legrendre_coeff(:,4) = 0.0625*[0 0 0 231 0 -315 0 105 0 -5];
    legrendre_coeff(:,5) = 0.0625*[0 0 429 0 -693 0 315 0 -35 0];
    legrendre_coeff(:,6) = 0.0078125*[0 6435 0 -12012 0 6930 0 -1260 0 35];
    legrendre_coeff(:,7) = 0.0078125*[12155 0 -25740 0 18018 0 -4620 0 315 0];

    derivative_coeff = [0 9 8 7 6 5 4 3 2 1]';
    legrendre_diff(2:end,:) = legrendre_coeff(1:end-1,:);
    legrendre_diff = legrendre_diff.*derivative_coeff;

    ind_map = [0 0 8 7 6 5 4 3 2];
    ind = ind_map(N);
    
    z = sort(roots(legrendre_diff(ind:end,N-2)));
    z = [-1; z; 1];
    
    % Evaluate Legendre polynomial at gll points
    pN_coeff = legrendre_coeff(:,N-2);
    pN_gll = pN_coeff(1)*z.^9 + pN_coeff(2)*z.^8 + pN_coeff(3)*z.^7 + ...
             pN_coeff(4)*z.^6 + pN_coeff(5)*z.^5 + pN_coeff(6)*z.^4 + ...
             pN_coeff(7)*z.^3 + pN_coeff(8)*z.^2 + pN_coeff(9)*z + ...
             pN_coeff(10);
    p_gll = pN_gll;
    pN_gll = pN_gll.^2;

    w = 2./((N+1)*N*pN_gll);
    
    Q = N + 1;
    D = zeros(Q,Q);

    for i = 1:Q
        for j = 1:Q
            if i == j && i == 1 && j == 1
                D(i,j) = -Q*(Q-1)/4;
            elseif i == Q && j == Q
                D(i,j) = Q*(Q-1)/4;
            elseif i >= 2 && i == j && j <= Q - 1
                D(i,j) = 0;
            elseif i >= 1 && j <= Q && i ~= j

                zij = [z(i) z(j)];

                ra = pN_coeff(1).*zij.^9 + pN_coeff(2).*zij.^8 + pN_coeff(3).*zij.^7 + ...
                         pN_coeff(4).*zij.^6 + pN_coeff(5).*zij.^5 + pN_coeff(6).*zij.^4 + ...
                         pN_coeff(7).*zij.^3 + pN_coeff(8).*zij.^2 + pN_coeff(9).*zij + ...
                         pN_coeff(10);

                D(i,j) = ra(1)/(ra(2)*(zij(1) - zij(2)));

            end
        end
    end 
    

end
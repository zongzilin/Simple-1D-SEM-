function x = std_element_mapping(E, N, Lx, z)
    % Each element length
    Le = Lx/E;
    Le_coord = Le:Le:Lx;
    
    % local to global mapping, do first element separatedly
    for i = 1:N+1
        x(i) = Le_coord(1)*(1+z(i))/2;
    end
    
    for iE = 2:E
        for iEl = 1:N+1
            x(iEl+(iE-1)*(N+1)) = Le_coord(iE-1)*(1-z(iEl))/2 + ...
                                Le_coord(iE)*(1+z(iEl))/2;
        end
    end
    x = x';

end
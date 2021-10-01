function L = spatial_derivatives(ws, ALPHA, BETA, GAMMA, M, DELTA_X)
    
    L = zeros(M, 1);
    
    if (GAMMA == 0)
        %% GAMMA == 0 case
        
        % Boundary near x = 0
        L(1) = 2 * (ws(2) - ws(1));
        
        % Bulk nodes
        L(2 : M - 1) = ws(1 : M - 2) - 2 * ws(2 : M - 1) + ws(3 : M);
        
        % Boundary near x = L
        L(M) = ws(M - 1) - 2 * ws(M);
        
        % Multiply by constant factor
        L = (BETA / (ALPHA * DELTA_X^2)) * L;
        
    else
        %% GAMMA > 0 case
        Cbeta = BETA * DELTA_X^2 / GAMMA;

        % Boundary near x = 0
        L(1) = -(6 + 2 * Cbeta) * ws(1) + 2 * (4 + Cbeta) * ws(2) - 2 * ws(3);
        L(2) = (4 + Cbeta) * ws(1) - (7 + 2 * Cbeta) * ws(2) + (4 + Cbeta) * ws(3) - ws(4);

        % Bulk nodes
        L(3 : M - 2) = ...
               - ws(1 : M - 4) ...
               + (4 + Cbeta) * ws(2 : M - 3) ...
               - (6 + 2 * Cbeta) * ws(3 : M - 2) ...
               + (4 + Cbeta) * ws(4 : M - 1) ...
               - ws(5 : M);

       % Boundary near x = L
       L(M - 1) = -ws(M - 3) + (4 + Cbeta) * ws(M - 2) - (6 + 2 * Cbeta) * ws(M - 1) + (4 + Cbeta) * ws(M);
       L(M) = -ws(M - 2) + (4 + Cbeta) * ws(M - 1) - (5 + 2 * Cbeta) * ws(M);

       % Multiply by constant factor
       L = (GAMMA / (ALPHA * DELTA_X^4)) * L;
    end

end
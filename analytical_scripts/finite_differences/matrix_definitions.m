function [A_mat, B_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T)


    if (GAMMA == 0) 
        Dbeta2 = (ALPHA * DELTA_X * DELTA_X) / (BETA * DELTA_T * DELTA_T);
        
        % A definition
        A_upper = -ones(M, 1);
        A_upper(2) = -2;
        
        A_main = (4 * Dbeta2 + 2) * ones(M, 1);
        
        A_lower = -ones(M, 1);
        
        A_mat = spdiags([A_lower A_main A_upper], -1:1, M, M);
        
        % B definition
        B_upper = 2 * ones(M, 1);
        B_upper(2) = 4;
        
        B_main = (8 * Dbeta2 - 4) * ones(M, 1);
        
        B_lower = 2 * ones(M, 1);
        
        B_mat = spdiags([B_lower B_main B_upper], -1:1, M, M);
        
    else
        Calpha = ALPHA * DELTA_X^4 / (GAMMA * DELTA_T^2);
        Cbeta = BETA * DELTA_X^2 / GAMMA;

        % A definition
        A_upper_upper = ones(M, 1);
        A_upper_upper(3) = 2;

        A_upper = (-Cbeta - 4) * ones(M, 1);
        A_upper(2) = -2 * Cbeta - 8;

        A_main = (4 * Calpha + 2 * Cbeta + 6) * ones(M, 1);
        A_main(2) = A_main(2) + 1;
        A_main(M) = A_main(M) - 1;

        A_lower = (-Cbeta - 4) * ones(M, 1);

        A_lower_lower = ones(M, 1);

        A_mat = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

        % B definition
        B_upper_upper = -2 * ones(M, 1);
        B_upper_upper(3) = 2 * B_upper_upper(3);

        B_upper = (2 * Cbeta + 8) * ones(M, 1);
        B_upper(2) = 2 * B_upper(2);
        B_main = (8 * Calpha - 4 * Cbeta - 12) * ones(M, 1);
        B_main(2) = B_main(2) - 2;
        B_main(M) = B_main(M) + 2;

        B_lower = (2 * Cbeta + 8) * ones(M, 1);

        B_lower_lower = -2 * ones(M, 1);

        B_mat = spdiags([B_lower_lower B_lower B_main B_upper B_upper_upper], -2:2, M, M);
    end
end
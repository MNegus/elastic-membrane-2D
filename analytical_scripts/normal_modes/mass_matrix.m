function [M, S] = mass_matrix(d, alpha, epsilon, L, N)

    %% Save lambdas
    lambdas = pi * (2 * (1 : N) - 1) / (2 * L);

    %% Calculate the S matrix
    % Bulk nodes for m != n
    S = (pi * d / epsilon) ...
        * (besselj(0, epsilon * d * lambdas)' * (lambdas .* besselj(1, epsilon * d * lambdas)) ...
        - (lambdas .* besselj(1, epsilon * d * lambdas))' * besselj(0, epsilon * d * lambdas)) ...
        ./ (lambdas.^2 - (lambdas.^2)');
    
    % Diagonals for m == n
    S(1 : N + 1 : end) = (pi * d^2 / 2) ...
        * (besselj(0, epsilon * d * lambdas).^2 + besselj(1, epsilon * d * lambdas).^2);
    
    %% Write the overall matrix
    M = alpha * eye(N) + S / L;
    
    %% Manual loop calculation for checking (keep commented)
%     Snn = @(n) (pi * d^2 / 2) * (besselj(0, epsilon * lambdas(n) * d)^2 + besselj(1, epsilon * lambdas(n) * d)^2);
%     Smn = @(m, n) (pi * d / epsilon) ...
%     * (lambdas(n) * besselj(0, epsilon * lambdas(m) * d) * besselj(1, epsilon * lambdas(n) * d) ...
%         - lambdas(m) * besselj(0, epsilon * lambdas(n) * d) * besselj(1, epsilon * lambdas(m) * d)) ...
%         / (lambdas(n)^2  - lambdas(m)^2);
%     A2 = zeros(N, N);
%     for m = 1 : N
%         for n = 1 : N
%             if m == n
%                A2(m, n) = Snn(n);
%             else
%                A2(m, n) = Smn(m, n);
%             end
%         end
%     end
%     max(max(S - A2))
end

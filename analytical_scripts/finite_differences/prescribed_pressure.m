function ps = prescribed_pressure(xs, alpha, epsilon, L, Ne)

    lambdas = pi * (2 * (1 : Ne) - 1) / (2 * L);
    p_ns = (alpha / epsilon^2) ./ sqrt(L * lambdas);
    ps = sum(p_ns .* cos(xs * lambdas), 2) / sqrt(L);
end
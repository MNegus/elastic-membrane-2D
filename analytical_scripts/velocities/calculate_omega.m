function [omega, Zs, Us] = calculate_omega()
    %% Find values of eta
    etas = jet_root_etas(0, 1024);

    %% Find the values of U(0, Z)
    Zs = zeros(size(etas));
    Us = zeros(size(etas));
    for m = 1 : length(etas)
        eta = etas(m);

        % Determine Z
        zeta = (1 / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
        Zs(m) = imag(zeta);

        % Save U
        Us(m) = real((1 + 1i * sqrt(eta)) / (1 - 1i * sqrt(eta)));
    end
    
    %% Find omega by integrating U^2
    omega = trapz(Zs, Us.^2);
end
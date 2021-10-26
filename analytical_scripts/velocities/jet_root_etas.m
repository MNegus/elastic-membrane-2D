function etas = jet_root_etas(no_points)

    %% Parameterises etas
    options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);
    eta0_guess = -0.01 + 1e-10;
    eta0 = fsolve(@(eta) (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1), eta0_guess, options);

    num_reals = 1e2;
    num_imags = 1e2;
    eta_reals = linspace(eta0, 1.1, num_reals);
    eta_imags = linspace(0, 0.5, num_imags);

    %% Finds x and z values
    xs = [];
    zs = [];
    eta_points = [];
    for eta_real = eta_reals
        for eta_imag = eta_imags
            eta = eta_real + 1i * eta_imag;
            zeta = (1 / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);

            xs(end + 1) = real(zeta);
            zs(end + 1) = imag(zeta);
            eta_points(end + 1) = eta;
        end 
    end

    %% Find query points in z
    zq = linspace(0, 1 + 4 / pi, no_points);
    xq = zeros(size(zq));

    etas = griddata(xs, zs, eta_points, xq, zq);

end
clear;

t = 0.1;
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 / sqrt(t);
J = @(t) (pi / 8) * d(t) / (d_t(t)^2);

%% Sets up figure
close(figure(1));
figure(1);
hold on;

%% Plots free surface
sigmas = 10.^linspace(-6, 1.1, 1e3);
surface_xs = (J(t) / pi) * (sigmas - log(sigmas) - 1);
surface_zs = J(t) * (1 + 4 * sqrt(sigmas) / pi);
plot(surface_xs, surface_zs);

%% Plot the sample points
options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);
eta0_guess = -0.01 + 1e-10;
eta0 = fsolve(@(eta) (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1), eta0_guess, options);
num_reals = 1e2;
num_imags = 1e2;
eta_reals = linspace(eta0, 2, num_reals);
eta_imags = linspace(0, 0.5, num_imags);
xs = [];
zs = [];
for eta_imag = eta_imags
    for eta_real = eta_reals
        eta = eta_real + 1i * eta_imag;
        zeta = (J(t) / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);

        xs(end + 1) = real(zeta);
        zs(end + 1) = imag(zeta);
    end 
    plot(surface_xs, surface_zs);
    hold on;
    scatter(xs, zs);
    hold off;
    drawnow;
    pause(0.01);
end

%% Find query points in z
no_points = 256;
etaq = jet_root_etas(no_points);
figure(1);
hold on;
xs_plot = zeros(size(etaq));
zs_plot = zeros(size(etaq));
for m = 1 : length(etaq)
    eta = etaq(m);
    zeta = (J(t) / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
        
    xs_plot(m) = real(zeta);
    zs_plot(m) = imag(zeta);
end
scatter(xs_plot, zs_plot, [], 'red');
xlim([-0.1, 0.1]);

%%



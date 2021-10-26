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



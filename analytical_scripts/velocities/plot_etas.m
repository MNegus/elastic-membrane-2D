clear;
close all;


%% Load in etas
X = 1;
etas = jet_root_etas(X, 1024);

%% Plot eta values
figure(1);

sigmas = real(etas);
xis = imag(etas);

zetas = (1 / pi) * (etas + 4 * 1i * sqrt(etas) - log(etas) + 1i * pi - 1);
xs = real(zetas);
zs = imag(zetas);

scatter(real(etas), imag(etas));

% Find the maximum point
[max_xi, idx] = max(xis);
max_sigma = sigmas(idx);
hold on;
scatter(max_sigma, max_xi);

% Find r and theta
zeta_max = max_sigma + 1i * max_xi;
arg = angle(zeta_max);
r = sqrt(max_sigma^2 + max_xi^2);
plot([0, r * cos(arg)], [0, r * sin(arg)]);
hold off;

thing = 1 + (2 / sqrt(r)) * sin(arg / 2) - (1 / r) * cos(arg)

% Find derivative
deriv = (xis(idx + 1) - xis(idx - 1)) / (sigmas(idx + 1) - sigmas(idx - 1))

%% (x, z) values


%% Plot z vs sigmas
figure(2);
plot(sigmas, zs);
hold on;
plot(sigmas, 4 * sqrt(sigmas - min(sigmas)));
hold off;

%% Plot derivative of zs with respect to sigma
figure(3);
plot(sigmas(1 : end - 1), diff(zs) ./ diff(sigmas));
hold on;
plot(sigmas, 1 ./ sqrt(sigmas - min(sigmas)));
hold off;



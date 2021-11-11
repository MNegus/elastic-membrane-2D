%% plot_energy.m

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();

% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIMESTEP = IMPACT_TIME / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

% Pressure type (composite or outer)
pressure_type = "composite";


%% Load in value of omega
omega_mat = load("omega.mat");
omega = omega_mat.omega;

%% Data dirs
analytical_parent_dir = "/media/michael/newarre/elastic_membrane/confirmation_data/gamma_varying_with_jet/analytical_data";

%% Stationary values
ds_stationary = 2 * sqrt(ts_analytical);
d_ts_stationary = 1 ./ sqrt(ts_analytical);
Js_stationary = pi * ds_stationary ./ (8 * d_ts_stationary.^2);
fluxes_stationary = 2 * Js_stationary .* d_ts_stationary.^3;
fluxes_stationary(1) = 0;
energy_stationary = cumtrapz(ts_analytical, fluxes_stationary);

%% Stationary jet energy comparison
close(figure(1));
figure(1);
hold on;

% Plot analytical value
plot(ts_analytical(2 : end), fluxes_stationary(2 : end));
MIN_CELL_SIZE = 6 / 2^12;
% Loads in Basilisk value
stationary_dir = "/media/michael/newarre/elastic_membrane/flux_with_pressure/fluxes";
flux_steps = 1290 : 10 : 3000;
basilisk_fluxes = zeros(size(flux_steps));
basilisk_flux_times = flux_steps * DELTA_T - IMPACT_TIME;
for idx = 1 : length(flux_steps)
	k = flux_steps(idx);
    flux_mat = dlmread(sprintf("%s/fluxes_%d.txt", stationary_dir, k));
    basilisk_fluxes(idx) = flux_mat(1, 2) * MIN_CELL_SIZE;
end

plot(basilisk_flux_times, basilisk_fluxes)





%% Jet energy comparison
close(figure(1));
figure(1);
hold on;

for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              analytical_parent_dir, ALPHA, BETA, GAMMA, pressure_type)
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
            
            % Extract d_ts and Js
            Js_mat = matfile(sprintf("%s/Js.mat", parameter_dir));
            Js = Js_mat.Js;
            
            d_ts_mat = matfile(sprintf("%s/d_ts.mat", parameter_dir));
            d_ts = d_ts_mat.d_ts;
            
            % Determine jet flux
            fluxes = (1 + omega) * Js .* d_ts.^3;
            
            % Determines jet energy
            jet_energy = EPSILON^2 * cumtrapz(ts_analytical, fluxes);
            
            % Plot line
            plot(ts_analytical, jet_energy, 'linewidth', 2, 'Displayname', displayname);
        end
    end
end
plot(ts_analytical, energy_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');
title("Analytical jet energy", "Interpreter", "latex", "Fontsize", 12);
legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
xlabel("$t$", 'Interpreter', 'Latex', 'Fontsize', 12);
ylabel("Total energy in jet", 'Interpreter', 'Latex', 'Fontsize', 12);


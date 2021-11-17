%% plot_energy.m
clear;

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();
no_params = length(ALPHAS);

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
analytical_parent_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping/alpha_varying";

%% Stationary values
ds_stationary = 2 * sqrt(ts_analytical);
d_ts_stationary = 1 ./ sqrt(ts_analytical);
Js_stationary = pi * ds_stationary ./ (8 * d_ts_stationary.^2);
fluxes_stationary = (1 + omega) * Js_stationary .* d_ts_stationary.^3;
fluxes_stationary(1) = 0;
energy_stationary = cumtrapz(ts_analytical, fluxes_stationary);

%% Jet energy comparison
close(figure(1));
figure(1);
hold on;

final_vals = zeros(no_params, 1);
for idx = 1 : no_params
    ALPHA = ALPHAS(idx);
    BETA = BETAS(idx);
    GAMMA = GAMMAS(idx);
    
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

    final_vals(idx) = jet_energy(end);
    
    % Plot line
    plot(ts_analytical, jet_energy, 'linewidth', 2, 'Displayname', displayname);
end
plot(ts_analytical, energy_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');
title("Analytical jet energy", "Interpreter", "latex", "Fontsize", 12);
legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
xlabel("$t$", 'Interpreter', 'Latex', 'Fontsize', 12);
ylabel("Total energy in jet", 'Interpreter', 'Latex', 'Fontsize', 12);

%%
close(figure(89));
figure(89);
plot(1 : no_params, final_vals, '-o');

%% Jet root height comparison
close(figure(2));
figure(2);
hold on;

for idx = 1 : no_params
    ALPHA = ALPHAS(idx);
    BETA = BETAS(idx);
    GAMMA = GAMMAS(idx);
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              analytical_parent_dir, ALPHA, BETA, GAMMA, pressure_type)
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
            
            % Extract d_ts and Js
            Js_mat = matfile(sprintf("%s/Js.mat", parameter_dir));
            Js = Js_mat.Js;
            
            
            % Plot line
            plot(ts_analytical(1 : end - 1), Js(1 : end - 1) * (1 + 4 / pi), 'linewidth', 2, 'Displayname', displayname);
%         end
end
plot(ts_analytical, Js_stationary * (1 + 4 / pi), 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
title("Analytical jet root height", "Interpreter", "latex", "Fontsize", 12);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
xlabel("$t$", 'Interpreter', 'Latex', 'Fontsize', 12);
ylabel("$(1 + 4 / \pi) J(t)$", 'Interpreter', 'Latex', 'Fontsize', 12);

%% Turnover point position
close(figure(3));
figure(3);
hold on;

for idx = 1 : no_params
    ALPHA = ALPHAS(idx);
    BETA = BETAS(idx);
    GAMMA = GAMMAS(idx);
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              analytical_parent_dir, ALPHA, BETA, GAMMA, pressure_type)
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
            
            % Extract d_ts and Js
            ds_mat = matfile(sprintf("%s/ds.mat", parameter_dir));
            ds = ds_mat.ds;
            
            
            % Plot line
            plot(ts_analytical(1 : end - 1), ds(1 : end - 1), 'linewidth', 2, 'Displayname', displayname);
end
plot(ts_analytical, ds_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
title("Analytical turnover point position", "Interpreter", "latex", "Fontsize", 12);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
xlabel("$t$", 'Interpreter', 'Latex', 'Fontsize', 12);
ylabel("$d(t)$", 'Interpreter', 'Latex', 'Fontsize', 12);

%% Turnover point velocity
close(figure(4));
figure(4);
hold on;

for idx = 1 : no_params
    ALPHA = ALPHAS(idx);
    BETA = BETAS(idx);
    GAMMA = GAMMAS(idx);
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              analytical_parent_dir, ALPHA, BETA, GAMMA, pressure_type)
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
            
            % Extract d_ts 
            d_ts_mat = matfile(sprintf("%s/d_ts.mat", parameter_dir));
            d_ts = d_ts_mat.d_ts;
            
            % Plot line
            plot(ts_analytical(1 : end - 1), d_ts(1 : end - 1), 'linewidth', 2, 'Displayname', displayname);

end
plot(ts_analytical, d_ts_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

legend("location", "northeast", "Interpreter", "latex", "Fontsize", 12);
title("Analytical turnover point velocity", "Interpreter", "latex", "Fontsize", 12);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
xlabel("$t$", 'Interpreter', 'Latex', 'Fontsize', 12);
ylabel("$d'(t) / \epsilon$", 'Interpreter', 'Latex', 'Fontsize', 12);
ylim([0, 20]);
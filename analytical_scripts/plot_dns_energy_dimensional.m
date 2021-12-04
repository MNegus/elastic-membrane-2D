%% plot_energy.m
clear;

%% Data definitions
master_dir = "/media/michael/newarre/elastic_membrane/basilisk_parameter_sweeping";

%% Dimensional parameters
rho = 998.0;
V = 5;
R = 1e-3;

energy_scale = rho * V^2 * R^2;
millimetre_scale = R / 1e-3;

[EPSILON, ~, ~, ~, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();

%% Derived parameters
% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIMESTEP = IMPACT_TIME / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_dimensional = 1000 * (R / V) * T_VALS;
NO_TIMESTEPS = length(T_VALS);
fluxes_steps = 1290 : 10 : 4000;

%% Loop over varying types
for varying = ["modulus"]
    %% Sets the parameters
    if varying == "thickness"
        ALPHAS = [1, 1.5, 2, 3, 4, 6, 8] / EPSILON^2;
        BETAS = zeros(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "tension"
        BETAS = [0, 10, 40, 160, 640, 2560, 10240] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "modulus"
        GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        ALPHAS = 1 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
    end
    % Parameters
    no_params = length(ALPHAS);

    %% Data dirs
    dns_parent_dir = sprintf("%s/%s_varying", master_dir, varying);

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
        parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/raw_data", ...
          dns_parent_dir, ALPHA, BETA, GAMMA);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
        
        %% DETERMINE THE FLUX, ONLY EVER NEEDS TO BE DONE ONCE
%         % Finds the fluxes
%         fluxes = zeros(NO_TIMESTEPS, 1);
%         for idx = fluxes_steps
%             idx
%             velocities_mat = readmatrix(sprintf("%s/velocities_%d-x_cell_0.txt", parameter_dir, idx));
%             zs = velocities_mat(:, 1);
%             fs = velocities_mat(:, 2);
%             us = velocities_mat(:, 3);
%             vs = velocities_mat(:, 4);
%             ps = velocities_mat(:, 5);
%             if (length(zs) <= 1)
%                 fluxes(idx) = 0;
%             else
%                 fluxes(idx) = trapz(zs, fs .* (0.5 * (us.^2 + vs.^2) + ps));
%             end
%         end
%         % Integrates the fluxes to get the energy
%         writematrix(non_dimensional_jet_energy, sprintf("%s/jet_energy.txt", parameter_dir));
%         jet_energy = energy_scale * non_dimensional_jet_energy;
        
        %% READ JET ENERGY IF SAVES BEFORE
        non_dimensional_jet_energy = readmatrix(sprintf("%s/jet_energy.txt", parameter_dir));
        jet_energy = energy_scale * non_dimensional_jet_energy;
        
        % Plot line
        plot(ts_dimensional, jet_energy, 'linewidth', 2, 'Displayname', displayname);
        drawnow;
    end
%     plot(ts_dimensional, energy_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');
    title("DNS jet energy", "Interpreter", "latex", "Fontsize", 12);
    legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    ylabel("Energy in jet (J/m)", 'Interpreter', 'Latex', 'Fontsize', 12);


    %% Jet root height comparison
    close(figure(2));
    figure(2);
    hold on;

    for idx = 1 : no_params
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);
        
        % Loads in parameters
        parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/raw_data", ...
          dns_parent_dir, ALPHA, BETA, GAMMA);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];

        % Extract ds
        turnover_points_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", parameter_dir));
        heights = turnover_points_mat(:, 3);

        % Plot line
        plot(ts_dimensional, millimetre_scale * heights, 'linewidth', 2, 'Displayname', displayname);
    end
%     plot(ts_dimensional, millimetre_scale * Js_stationary * (1 + 4 / pi), 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

    legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
    title("DNS jet root height", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    ylabel("Jet root height (mm)", 'Interpreter', 'Latex', 'Fontsize', 12);

    %% Turnover point position
    close(figure(3));
    figure(3);
    hold on;

    for idx = 1 : no_params
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);
        
        % Loads in parameters
        parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/raw_data", ...
          dns_parent_dir, ALPHA, BETA, GAMMA);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];

        % Extract ds
        turnover_points_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", parameter_dir));
        ds = turnover_points_mat(:, 2);

        % Plot line
        plot(ts_dimensional, millimetre_scale * ds, 'linewidth', 2, 'Displayname', displayname);
    end
%     plot(ts_dimensional, millimetre_scale * ds_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

    legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
    title("DNS turnover point position", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    ylabel("Jet root $x^*$ position (mm)", 'Interpreter', 'Latex', 'Fontsize', 12);

    %% Turnover point velocity
    close(figure(4));
    figure(4);
    hold on;

    for idx = 1 : no_params
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);
        
        % Loads in parameters
        parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/raw_data", ...
          dns_parent_dir, ALPHA, BETA, GAMMA);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];

        % Extract ds
        turnover_points_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", parameter_dir));
        d_ts = turnover_points_mat(:, 4);

        % Plot line
        plot(ts_dimensional(2 : end - 1), V * d_ts(2 : end - 1), 'linewidth', 2, 'Displayname', displayname);

    end
%     plot(ts_dimensional(1 : end - 1), V * d_ts_stationary(1 : end - 1), 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

    legend("location", "northeast", "Interpreter", "latex", "Fontsize", 12);
    title("Analytical turnover point velocity", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);

    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    ylabel("Jet root $x$ velocity (m/s)", 'Interpreter', 'Latex', 'Fontsize', 12);
    ylim([0, 20]);
end
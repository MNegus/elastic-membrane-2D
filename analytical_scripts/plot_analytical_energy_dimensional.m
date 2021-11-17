%% plot_energy.m
clear;
close all;

%% Data definitions
master_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping";

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
ts_analytical = (0 : DELTA_T : T_MAX - IMPACT_TIME);
ts_dimensional = 1000 * (R / V) * ts_analytical;

% Pressure type (composite or outer)
pressure_type = "composite";

% Load in omega
omega_mat = load("omega.mat");
omega = omega_mat.omega;

%% Colour mapping
no_params = 7;
colormap("jet")
cmap = colormap;
color_idxs = floor(linspace(1, length(cmap), no_params));
    
%% Stationary values
ds_stationary = 2 * sqrt(ts_analytical);
d_ts_stationary = 1 ./ sqrt(ts_analytical);
Js_stationary = pi * ds_stationary ./ (8 * d_ts_stationary.^2);
fluxes_stationary = (1 + omega) * Js_stationary .* d_ts_stationary.^3;
fluxes_stationary(1) = 0;
energy_stationary = energy_scale * cumtrapz(ts_analytical, fluxes_stationary);

%% Loop over varying types
for varying = ["alpha", "alpha", "beta", "gamma"]
    %% Sets the parameters
    if varying == "alpha"
        ALPHAS = [1, 1.5, 2, 3, 4, 6, 8] / EPSILON^2;
        BETAS = zeros(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "beta"
        BETAS = [0, 10, 40, 160, 640, 2560, 10240] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "gamma"
        GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        ALPHAS = 2 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
    end
    % Parameters
    no_params = length(ALPHAS);

    %% Data dirs
    analytical_parent_dir = sprintf("%s/%s_varying", master_dir, varying);

    %% Set up figure
    layout = tiledlayout(2, 4);
    x0=400;
    y0=400;
    width=1200;
    height=601;
    set(gcf,'position',[x0,y0,width,height])
    
    %% Turnover point position
    nexttile(1, [2, 1]);
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
        plot(ts_dimensional(1 : end - 1), millimetre_scale * ds(1 : end - 1), 'linewidth', 2, ...
            'color', cmap(color_idxs(idx), :));
    end
    plot(ts_dimensional, millimetre_scale * ds_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

%     legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
%     title("Analytical turnover point position", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    xlim([0, 1.1 * max(ts_dimensional)]);
    ylabel("Jet root position (mm)", 'Interpreter', 'Latex', 'Fontsize', 12);
    
    %% Jet root height comparison
    nexttile(2, [2, 1]);
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
        plot(ts_dimensional(1 : end - 1), millimetre_scale * Js(1 : end - 1) * (1 + 4 / pi), 'linewidth', 2, ...
            'color', cmap(color_idxs(idx), :));
    end
    plot(ts_dimensional, millimetre_scale * Js_stationary * (1 + 4 / pi), 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

%     legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
%     title("Analytical jet root height", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    xlim([0, 1.1 * max(ts_dimensional)]);
    ylabel("Jet root height (mm)", 'Interpreter', 'Latex', 'Fontsize', 12);
    
    %% Turnover point velocity
    nexttile(3, [2, 1]);
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
        plot(ts_dimensional(2 : end - 1), V * d_ts(2 : end - 1), 'linewidth', 2, ...
            'color', cmap(color_idxs(idx), :));

    end
    plot(ts_dimensional(1 : end - 1), V * d_ts_stationary(1 : end - 1), 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');

%     legend("location", "northeast", "Interpreter", "latex", "Fontsize", 12);
%     title("Analytical turnover point velocity", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);

    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    xlim([0, 1.1 * max(ts_dimensional)]);
    ylabel("Jet root velocity (m/s)", 'Interpreter', 'Latex', 'Fontsize', 12);
    ylim([6, 20]);
    
    %% Jet energy comparison
    nexttile(4, [2, 1]);
    hold on;

    final_vals = zeros(no_params, 1);
    for idx = 1 : no_params
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);

        % Loads in parameters
        parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
          analytical_parent_dir, ALPHA, BETA, GAMMA, pressure_type);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];

        % Extract d_ts and Js
        Js_mat = matfile(sprintf("%s/Js.mat", parameter_dir));
        Js = Js_mat.Js;

        d_ts_mat = matfile(sprintf("%s/d_ts.mat", parameter_dir));
        d_ts = d_ts_mat.d_ts;

        % Determine jet flux
        fluxes = (1 + omega) * Js .* d_ts.^3;

        % Determines jet energy
        jet_energy = energy_scale * EPSILON^2 * cumtrapz(ts_analytical, fluxes);

        final_vals(idx) = jet_energy(end);

        % Plot line
        plot(ts_dimensional, jet_energy, 'linewidth', 2, ...
            'color', cmap(color_idxs(idx), :));
    end
    plot(ts_dimensional, energy_stationary, 'linewidth', 2, 'Displayname', "Stationary", 'color', 0.5 * [1 1 1], 'linestyle', '--');
%     title("Analytical jet energy", "Interpreter", "latex", "Fontsize", 12);
%     legend("location", "northwest", "Interpreter", "latex", "Fontsize", 12);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 12);
    xlabel("$t^*$ (ms)", 'Interpreter', 'Latex', 'Fontsize', 12);
    xlim([0, 1.1 * max(ts_dimensional)]);
    ylabel("Energy in jet (J/m)", 'Interpreter', 'Latex', 'Fontsize', 12);


    
    set(gcf,'position',[x0,y0,width,height])
    drawnow;
    exportgraphics(gcf, sprintf("%s_varying_quantities.png", varying));
    pause(1);
    
end



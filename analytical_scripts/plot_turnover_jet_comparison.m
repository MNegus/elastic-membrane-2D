%% plot_maximum_deflection.m
clear;
close all;

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE, IMPACT_TIME] ...
    = parameters();

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
timesteps = 1 : 2000;
ts_analytical = ts_analytical(timesteps);

fontsize = 22;

% Timestep to plot
timestep = ceil((IMPACT_TIME + 0.2) / DELTA_T) + 1;

%% Data dirs
master_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping";

%%
tiledlayout(2, 3, 'TileSpacing','Compact', 'Padding', 'Compact');

varying_strs = ["alpha", "beta", "gamma"]
for var_idx = [1, 2, 3]
    varying = varying_strs(var_idx);
    
    parent_dir = sprintf("%s/%s_varying", master_dir, varying);
    
    %% Sets the parameters
    if varying == "alpha"
        ALPHAS = 2.^[0, 0.5, 1, 1.5, 2, 2.5, 3] / EPSILON^2;
        BETAS = ones(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * ones(size(ALPHAS)) * EPSILON^2;
        independent_param = ALPHAS;
    elseif varying == "beta"
        BETAS = [5, 10, 20, 40, 80, 160, 320, 640, 1280] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
        independent_param = BETAS;
    elseif varying == "gamma"
        GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        ALPHAS = 1 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
        independent_param = GAMMAS;
    end
    
    no_params = length(ALPHAS);
    color_mags = linspace(0.75, 0, no_params + 1);
    
    %% Plots turnover points
    nexttile(var_idx);
    hold on;
    grid on;
    
    min_final = 1e6;
    for idx = 1 : no_params
        % Determine parameters
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);

        param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);

        % Load in ds
        ds_mat = matfile(sprintf("%s/finite_differences/composite/ds.mat", param_dir));
        ds = ds_mat.ds;

        % Plot ds
        plot(ts_analytical, ds(timesteps), ...
            'color', color_mags(idx) * [1 1 1], 'linewidth', 1.5);
        
        % Finds min final
        min_final = min(min_final, ds(timesteps(end)));
    end
    % Plot stationary value
    ds_stationary = 2 * sqrt(ts_analytical);
    plot(ts_analytical, ds_stationary, ...
        'color', 'black', 'linewidth', 2, 'linestyle', '--');
    
    % Figure properties
    xmax = 0.25;
    xlim([0, xmax]);
    set(gca,'XTickLabel',[]);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    
    if varying == "alpha"
        title("$\beta = 0$, $\gamma = 2$", "interpreter", "latex", "fontsize", fontsize);
        ylabel("$d_0(t)$", "interpreter", "latex", "fontsize", fontsize);
        annotation('arrow',[0.335 0.335], [0.8566 0.9112]);
        text(0.222, 0.87, "$\alpha$", "interpreter", "latex", "fontsize", fontsize);
    elseif varying == "beta"
        title("$\alpha = 1$, $\gamma = 2$", "interpreter", "latex", "fontsize", fontsize);
        set(gca,'YTickLabel',[]);
        annotation('arrow',[0.626 0.626], [0.8566 0.9112]);
        text(0.222, 0.87, "$\beta$", "interpreter", "latex", "fontsize", fontsize);
    else
        title("$\alpha = 1$, $\beta = 0$", "interpreter", "latex", "fontsize", fontsize);
        set(gca,'YTickLabel',[]);
        annotation('arrow',[0.915 0.915], [0.8566 0.9112]);
        text(0.222, 0.87, "$\gamma$", "interpreter", "latex", "fontsize", fontsize);
    end
    
    %% Plots jet root height
    nexttile(var_idx + 3);
    hold on;
    grid on;
    
    for idx = 1 : no_params
        % Determine parameters
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);

        param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);

        % Load in ds
        Js_mat = matfile(sprintf("%s/finite_differences/composite/Js.mat", param_dir));
        Js = Js_mat.Js;

        % Plot ds
        plot(ts_analytical, Js(timesteps), ...
            'color', color_mags(idx) * [1 1 1], 'linewidth', 1.5);
    end

    % Plot stationary value 
    plot(ts_analytical, pi * ts_analytical.^(3/2) / 4, ...
        'color', 'black', 'linewidth', 2, 'linestyle', '--');
    
    % Figure properties
    xlim([0, xmax]);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    xlabel("$t$", "interpreter", "latex", "Fontsize", fontsize);
    
    if (var_idx) > 1
        set(gca,'YTickLabel',[]);
    end
    
    if varying == "alpha"
        ylabel("$J(t)$", "interpreter", "latex", "fontsize", fontsize);
        annotation('arrow',[0.335 0.335], [0.287 0.428]);
        text(0.222, 0.0665, "$\alpha$", "interpreter", "latex", "fontsize", fontsize);
    elseif varying == "beta"
        annotation('arrow',[0.626 0.626], [0.287 0.428]);
        text(0.222, 0.0665, "$\beta$", "interpreter", "latex", "fontsize", fontsize);
    else
        annotation('arrow',[0.915 0.915], [0.287 0.428]);
        text(0.222, 0.0665, "$\gamma$", "interpreter", "latex", "fontsize", fontsize);
    end
end

%% Figure settings
x0=400;
y0=400;
width=800;
height=920;

set(gcf,'position',[x0,y0,width,height])

exportgraphics(gcf, "figures/turnover_jet_comparison.png", "Resolution", 300);
savefig(gcf, "figures/turnover_jet_comparison.fig");
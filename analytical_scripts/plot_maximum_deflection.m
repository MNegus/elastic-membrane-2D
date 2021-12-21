%% plot_maximum_deflection.m
clear;
close all;

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

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

fontsize = 22;

% Timestep to plot
timestep = ceil((IMPACT_TIME + 0.2) / DELTA_T) + 1;

%% Data dirs
master_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping";

%%
tiledlayout(1, 3, 'TileSpacing','Compact', 'Padding', 'Compact');

for varying = ["alpha", "beta", "gamma"]
    nexttile;
    hold on;
    parent_dir = sprintf("%s/%s_varying", master_dir, varying);
    
    %% Sets the parameters
    if varying == "alpha"
        ALPHAS = 2.^[0, 0.5, 1, 1.5, 2, 2.5, 3] / EPSILON^2;
        BETAS = zeros(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * ones(size(ALPHAS)) * EPSILON^2;
        independent_param = ALPHAS;
    elseif varying == "beta"
        BETAS = [5, 10, 20, 40, 80, 160, 320, 640, 1280] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2  * ALPHAS).^3 * EPSILON^2;
        independent_param = BETAS;
    elseif varying == "gamma"
%         GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        GAMMAS = [0.5, 1, 2, 4, 8, 16, 32] * EPSILON^2;
        ALPHAS = 2 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
        independent_param = GAMMAS;
    end
    
    no_params = length(ALPHAS);
    colors = ones(no_params, 3);
    
    color_idxs = floor(linspace(1, length(cmap), no_params));
    
    for q = 1 : no_params
        colors(q, :) = cmap(color_idxs(q), :);
    end
    
    %% Plots max deflection
    deflections = zeros(no_params, 1);
    for idx = 1 : no_params
       ALPHA = ALPHAS(idx);
       BETA = BETAS(idx);
       GAMMA = GAMMAS(idx);

       param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);

       % Loads in w values
       ws_mat = matfile(sprintf("%s/finite_differences/composite/w_%d.mat", param_dir, timestep - IMPACT_TIMESTEP));
       ws = EPSILON^2 * ws_mat.w_next;
       idx
       % Saves deflection at x = 0
       deflections(idx) = ws(1);
    end

    %% Plot deflections as a function of ALPHA
    plot(independent_param, deflections, 'color', 'black');
    sz = 75;
    scatter(independent_param, deflections, sz, colors, 'filled');
    set(gca, 'xscale', 'log');
%     set(gca, 'yscale', 'log');
    xticks(independent_param);
    ylim([0.01, 0.06]);
    
    
    %% Labels
    if varying == "alpha"
        title("$\beta = 0$, $\gamma = 2$", "interpreter", "latex", "fontsize", fontsize);
        xlim([2^(-0.5), 2^3.5]);
        xticks([1, 2, 4, 8]);
        xlabel("$\alpha$", "interpreter", "latex", "fontsize", fontsize);
        ylabel("max($w(x, t)$)", "interpreter", "latex", "fontsize", fontsize);
    elseif varying == "beta"
        title("$\alpha = 1$, $\gamma = 2$", "interpreter", "latex", "fontsize", fontsize);
        xticks(10 * [0.5, 2, 8, 32, 128]);
        xlim([10 * 2^(-2), 10 * 2^8]);
        set(gca,'YTickLabel',[]);
        xlabel("$\beta$", "interpreter", "latex", "fontsize", fontsize);
    else
        title("$\alpha = 2$, $\beta = 0$", "interpreter", "latex", "fontsize", fontsize);
%         xticks(10.^[-5, -4, -3, -2, -1, 0, 1, 2]);
%         xticks(10.^(-4 : 2 : 2));
        xlim([0.25, 64]);
        set(gca,'YTickLabel',[]);
        xlabel("$\gamma$", "interpreter", "latex", "fontsize", fontsize);
    end
    
    %% Figure properties
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize, 'fontsize', fontsize);
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
end

%% Figure settings
x0=400;
y0=400;
width=1200;
height=400;

set(gcf,'position',[x0,y0,width,height])

exportgraphics(gcf, "figures/deflection_comparison.png", "Resolution", 300);
savefig(gcf, "figures/deflection_comparison.fig");




%% plot_multiple_fd.m

close all;
clear;

%% Dimensional parameters
R = 1e-3;
V = 5;


%% Parameters
[EPSILON, ~, ~, ~, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

no_params = 7;

% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0;
T_VALS = 0 : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

% Pressure type (composite or outer)
pressure_type = "composite";


%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping";

%% Sets up figure
layout = tiledlayout(3, 3);
% layout.Padding = 'compact';
fontsize = 20;

x0=400;
y0=400;
width=1200;
height=801;

set(gcf,'position',[x0,y0,width,height])
%% Set up animated lines
nexttile(1, [1, 3]);
alpha_varying_lines = arrayfun(@(x) animatedline(), 1:no_params);
xlim([0, 6]);
xticks([0 : 0.5 : 6]);
ylim([-0.08, 0.02]);
ylabel("$-w^*(x^*, t^*)$ (mm)", "interpreter", "latex", "Fontsize", fontsize);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
set(gca,'Xticklabel',[]) 
grid on;

nexttile(4, [1, 3]);
beta_varying_lines = arrayfun(@(x) animatedline(), 1:no_params);
xlim([0, 6]);
xticks([0 : 0.5 : 6]);
ylim([-0.08, 0.02]);
ylabel("$-w^*(x^*, t^*)$ (mm)", "interpreter", "latex", "Fontsize", fontsize);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
set(gca,'Xticklabel',[]) 
grid on;

nexttile(7, [1, 3]);
gamma_varying_lines = arrayfun(@(x) animatedline(), 1:no_params);
xlim([0, 6]);
xticks([0 : 0.5 : 6]);
ylim([-0.08, 0.02]);
xlabel("$x^*$ (mm)", "interpreter", "latex", "Fontsize", fontsize);
ylabel("$-w^*(x^*, t^*)$ (mm)", "interpreter", "latex", "Fontsize", fontsize);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
grid on;

colormap("jet")
cmap = colormap;
color_idxs = floor(linspace(1, length(cmap), no_params));
for idx = 1 : no_params
    for line = [alpha_varying_lines(idx), beta_varying_lines(idx), gamma_varying_lines(idx)]
        color_idxs(idx);
        line.Color = cmap(color_idxs(idx), :);
        line.LineWidth = 2;
    end
end

%% Writer object setup
membrane_obj = VideoWriter('membrane.avi');
membrane_obj.FrameRate = 15;
open(membrane_obj);

%% Loops over time
for k = 1 : 20 : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t
    
    if (t <= 0) 
        continue 
    end
        
    varying_idx = 0;
    for varying = ["alpha", "beta", "gamma"]
        varying_idx = varying_idx + 1;
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
        
        
        
        %% Loops over parameters and plots
        varying_idx;
        nexttile(3 * (varying_idx - 1) + 1, [1, 3]);
        
        for idx = 1 : no_params
            ALPHA = ALPHAS(idx);
            BETA = BETAS(idx);
            GAMMA = GAMMAS(idx);

            % Loads in parameters
            parameter_dir = sprintf("%s/%s_varying/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              parent_dir, varying, ALPHA, BETA, GAMMA, pressure_type);
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];


            ws_mat = matfile(sprintf("%s/w_%d.mat", parameter_dir, k - IMPACT_TIMESTEP));
            ws = EPSILON^2 * ws_mat.w_next;

            % w plot   
            if varying == "alpha"
                title(sprintf("$t^*$ = %.4f (ms)", t / V), "Interpreter", "latex"); 
                
                clearpoints(alpha_varying_lines(idx));
                addpoints(alpha_varying_lines(idx), xs, -ws);
            elseif varying == "beta"
                clearpoints(beta_varying_lines(idx));
                addpoints(beta_varying_lines(idx), xs, -ws);
            elseif varying == "gamma"
                clearpoints(gamma_varying_lines(idx));
                addpoints(gamma_varying_lines(idx), xs, -ws);
            end
        end
    end
    
    f1 = getframe(gcf);
    writeVideo(membrane_obj, f1);
    
    drawnow;
end
close(membrane_obj);
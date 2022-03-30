%% plot_multiple_sna.m

close all;
clear;

%% Parameters
[~, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();
ALPHAS = 2.004;
BETAS = 0;
GAMMAS = 1069;

no_params = length(ALPHAS);


% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;


%% Data dirs
% parent_dir = "/media/michael/newarre/elastic_membrane/basilisk_parameter_sweeping/modulus_varying";
parent_dir = "/scratch/negus/realisticParamTests";



%% Loops over time
% close all;
close(figure(1));
figure(1);
for k = 1 : 100 : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t
    
    %% Loops over parameters and plots 
    % Resets holds
    subplot(3, 1, 1);
    hold off
    
    subplot(3, 1, 2);
    hold off;
    
    subplot(3, 1, 3);
    hold off;
    

    for idx = 1 : no_params
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);
                
        % Loads in parameters
%         parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/raw_data", ...
%           parent_dir, ALPHA, BETA, GAMMA);
        parameter_dir = sprintf("%s/ALPHA-%g_BETA-%g_GAMMA-%g/raw_data", ...
                  parent_dir, ALPHA, BETA, GAMMA);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];

        % w plot
        ws_mat = readmatrix(sprintf("%s/w_%d.txt", parameter_dir, k - IMPACT_TIMESTEP));
        [xs_sort, sort_idxs] = sort(ws_mat(:, 1));
        ws = ws_mat(sort_idxs, 2);
        
        subplot(3, 1, 1);
        plot(xs_sort, ws, 'Linewidth', 2, 'Displayname', displayname);
        hold on;
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

        % w_t plot
        w_ts_mat = readmatrix(sprintf("%s/w_deriv_%d.txt", parameter_dir, k - IMPACT_TIMESTEP));
        [xs_sort, sort_idxs] = sort(w_ts_mat(:, 1));
        w_ts = w_ts_mat(sort_idxs, 2);
        
        subplot(3, 1, 2);
        plot(xs_sort, w_ts, 'Linewidth', 2, 'Displayname', displayname);
        hold on;
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

        % p plot
        ps_mat = readmatrix(sprintf("%s/p_%d.txt", parameter_dir, k - IMPACT_TIMESTEP));
        [xs_sort, sort_idxs] = sort(ps_mat(:, 1));
        ps = ps_mat(sort_idxs, 2);
        
        subplot(3, 1, 3);
        plot(xs_sort, ps, 'Linewidth', 2, 'Displayname', displayname);
        hold on;
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$p(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
    end
    
    
    %% Figure settings
    x0=400;
    y0=400;
    width=1200;
    height=800;

    set(gcf,'position',[x0,y0,width,height])
    drawnow;

end
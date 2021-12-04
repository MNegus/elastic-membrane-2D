%% plot_multiple_fd.m

close all;
clear;

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

no_params = length(ALPHAS);


% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0;
T_VALS = 0 : DELTA_T : T_MAX;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

% Pressure type (composite or outer)
pressure_type = "composite";


%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping/gamma_varying";

%% Turnover point comparison
close(figure(1));
figure(1);
hold on;

for idx = 1 : no_params
    ALPHA = ALPHAS(idx);
    BETA = BETAS(idx);
    GAMMA = GAMMAS(idx);
    
    % Loads in parameters
    parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
      parent_dir, ALPHA, BETA, GAMMA, pressure_type)
    displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];

    % FD turnover points
    fd_comp_mat = matfile(sprintf("%s/ds.mat", parameter_dir));
    ds_comp = fd_comp_mat.ds;

    % Plot line
    plot(ts_analytical(1 : end - 1), ds_comp(1 : end - 1), 'linewidth', 2, 'Displayname', displayname);

end
legend("location", "northeast");

%% Loops over time
% close all;
close(figure(2));
figure(2);
for k = 1 : 100 : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t
    
    if (t <= 0) 0.125 / DELTA_T
        continue 
    end
        
    
    %% Loops over parameters and plots 
    % Resets holds
    subplot(3, 1, 1);
    hold off
    
    subplot(3, 1, 2);
    hold off;
    
    subplot(3, 1, 3);
    hold off;
    
%     figure(1);
%     hold off;

    for idx = 1 : no_params
        ALPHA = ALPHAS(idx);
        BETA = BETAS(idx);
        GAMMA = GAMMAS(idx);
                
        % Loads in parameters
        parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
          parent_dir, ALPHA, BETA, GAMMA, pressure_type);
        displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];


        ws_mat = matfile(sprintf("%s/w_%d.mat", parameter_dir, k - IMPACT_TIMESTEP));
        ws = EPSILON^2 * ws_mat.w_next;

        w_ts_mat = matfile(sprintf("%s/w_t_%d.mat", parameter_dir, k - IMPACT_TIMESTEP));
        w_ts = EPSILON^2 * w_ts_mat.w_t;

        ps_mat = matfile(sprintf("%s/p_%d.mat", parameter_dir, k - IMPACT_TIMESTEP));
        ps = ps_mat.p;

        % w plot
        subplot(3, 1, 1);
        plot(xs, ws, 'Linewidth', 2, 'Displayname', displayname);
        hold on;

%                 xlim([0, 2]);
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%                 legend("interpreter", "latex");
        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

        % w_t plot
        subplot(3, 1, 2);
        plot(xs, w_ts, 'Linewidth', 2, 'Displayname', displayname);
        hold on;

%                 xlim([0, 2]);
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%                 legend("interpreter", "latex");
        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

        % p plot
        subplot(3, 1, 3);
        plot(xs, ps, 'Linewidth', 2, 'Displayname', displayname);
        hold on;

%                 xlim([0, 2]);
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$p(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%                 legend("interpreter", "latex");
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
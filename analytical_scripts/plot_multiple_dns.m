%% plot_multiple_fd.m



%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;


%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/confirmation_data/gamma_varying/basilisk_data";

%% Loops over time
% close all;
for k = IMPACT_TIMESTEP : 100 : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t
    
    if (t <= 0) 
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
    
    for ALPHA = ALPHAS
        for BETA = BETAS
            for GAMMA = GAMMAS
                
                % Loads in parameters
                parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", ...
                  parent_dir, ALPHA, BETA, GAMMA);
                displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
                
                % Load in ws
                w_membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", parameter_dir, k - 1));
                unsorted_xs = w_membrane_mat(:, 1);
                unsorted_ws = w_membrane_mat(:, 2);
                [sorted_xs, idxs] = sort(unsorted_xs);
                ws = unsorted_ws(idxs);
                
                % Load in w_ts
                w_t_membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", parameter_dir, k - 1));
                unsorted_xs = w_t_membrane_mat(:, 1);
                unsorted_w_ts = w_t_membrane_mat(:, 2);
                [sorted_xs, idxs] = sort(unsorted_xs);
                w_ts = unsorted_w_ts(idxs);
                
                % Load in ps
                p_membrane_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", parameter_dir, k - 1));
                unsorted_xs = p_membrane_mat(:, 1);
                unsorted_ps = p_membrane_mat(:, 2);
                [sorted_xs, idxs] = sort(unsorted_xs);
                ps = unsorted_ps(idxs);
                
                % w plot
                subplot(3, 1, 1);
                plot(sorted_xs, ws, 'Linewidth', 2, 'Displayname', displayname);
                hold on;
                
%                 xlim([0, 2]);
                xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
                ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
                set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%                 legend("interpreter", "latex");
                title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
                
                % w_t plot
                subplot(3, 1, 2);
                plot(sorted_xs, w_ts, 'Linewidth', 2, 'Displayname', displayname);
                hold on;
                
%                 xlim([0, 2]);
                xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
                ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
                set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%                 legend("interpreter", "latex");
                title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
                
                % p plot
                subplot(3, 1, 3);
                plot(sorted_xs, ps, 'Linewidth', 2, 'Displayname', displayname);
                hold on;
                
%                 xlim([0, 2]);
                xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
                ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
                set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%                 legend("interpreter", "latex");
                title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
                
            end
        end
    end
    
    
    %% Figure settings
    x0=400;
    y0=400;
    width=1200;
    height=800;

    set(gcf,'position',[x0,y0,width,height])
    drawnow;
    
    
    
    
    
end
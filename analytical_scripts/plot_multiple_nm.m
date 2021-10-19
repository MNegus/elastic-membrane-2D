%% plot_multiple_nm.m

%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0;
IMPACT_TIMESTEP = 1;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;


%% Data dirs
parent_dir = "/home/negus/Desktop/alpha_vary_test";


%% Turnover point comparison
close(figure(1));
figure(1);
hold on;
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/normal_modes", ...
              parent_dir, ALPHA, BETA, GAMMA);
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
          
            % Normal modes turnover points
            nm_mat = matfile(sprintf("%s/ds.mat", parameter_dir));
            ds_nm = nm_mat.ds;
            
            % Plot line
            plot(ts_analytical, ds_nm, 'linewidth', 2, 'Displayname', displayname);
            pause(1);
        end
    end
end
legend("location", "northeast");

%% Loops over time
% close all;
for k = 1 : 100 : length(T_VALS)
% for k = IMPACT_TIMESTEP + [51, 501, 1001, 2001] 
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
                
                N_mat = matfile(sprintf("%s/normal_modes/N.mat", parameter_dir));
                N = N_mat.N

                as_mat = matfile(sprintf("%s/normal_modes/as.mat", parameter_dir));
                as = as_mat.as;

                a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", parameter_dir));
                a_ts = a_ts_mat.a_ts;

                q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", parameter_dir));
                q_ts = q_ts_mat.q_ts;
                
                nm_mat = matfile(sprintf("%s/normal_modes/ds.mat", parameter_dir));
                ds_nm = nm_mat.ds;
                
                % Normal modes
                [ws_nm, w_ts_nm, ps_nm] ...
                    = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
                    a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
                    ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON); 
            
                % w plot
                subplot(3, 1, 1);
                plot(xs, ws_nm, 'Linewidth', 2, 'Displayname', displayname);
                hold on;
                
%                 xlim([0, 2]);
                xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
                ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
                set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
                legend("interpreter", "latex");
                title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
                
                % w_t plot
                subplot(3, 1, 2);
                plot(xs, w_ts_nm, 'Linewidth', 2, 'Displayname', displayname);
                hold on;
                
%                 xlim([0, 2]);
                xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
                ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
                set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
                legend("interpreter", "latex");
%                 title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
                
                % p plot
                subplot(3, 1, 3);
                plot(xs, ps_nm, 'Linewidth', 2, 'Displayname', displayname);
                hold on;
                
%                 xlim([0, 2]);
                xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
                ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
                set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
                legend("interpreter", "latex");
%                 title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
                
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
    
    exportgraphics(gcf, sprintf("overlayed_membrane_%d.png", k), 'resolution', 300);
    
    
    
    
    
end
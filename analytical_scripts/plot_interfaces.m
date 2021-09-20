%% plot_interfaces.m

addpath("interface_analysis");

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


%% Turnover point comparison
close(figure(1));
figure(1);
hold on;

close(figure(2));
figure(2);
hold on;
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", ...
              parent_dir, ALPHA, BETA, GAMMA)
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
          
            % DNS turnover points
            dns_mat = dlmread(sprintf("%s/turnover_points.txt", parameter_dir));
            ds_dns = dns_mat(:, 2);
            ts_dns = (dns_mat(:, 1) - IMPACT_TIMESTEP) * DELTA_T;

            % Find a polynomial representation of ds for positive time
            idx = find(dns_mat(:, 1) == IMPACT_TIMESTEP);
            ds_pos = ds_dns(idx : end);
            ds_dns_sq = ds_pos.^2;
            p = polyfit(ts_analytical, ds_dns_sq, 100);
            ds_dns_poly = sqrt(polyval(p, ts_analytical)); 

            % Differentiate the polynomial to get d'(t)
            p_dir = polyder(p);
            d_ts_dns = (1 ./ (2 *  ds_pos)) .* polyval(p_dir, ts_analytical)';
            
            %% Plots d
            figure(1);
            plot(ts_dns, ds_dns.^2, 'linewidth', 2, 'displayname', displayname);
%             pause(1);
            
%             %% Plots d_ts
%             figure(2);
%             plot(ts_analytical, d_ts_dns, 'displayname', displayname);
            
        end
    end
end
xlabel("$t$", "interpreter", "latex", "Fontsize", 18);
ylabel("$d(t)^2$", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
legend("interpreter", "latex", "location", "northwest");
grid on;
x0=400;
y0=400;
width=1200;
height=500;

    set(gcf,'position',[x0,y0,width,height])
    drawnow;


%% Loops over time
% close all;
for k = 2000
    %% Updates time
    t = T_VALS(k);
    t
    
    if (t <= 0) 
        continue 
    end
        
    
    %% Loops over parameters and plots 
    % Resets holds
    figure(1);
    hold off;
    
    
    for ALPHA = ALPHAS
        for BETA = BETAS
            for GAMMA = GAMMAS
                
                % Loads in parameters
                parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", ...
                  parent_dir, ALPHA, BETA, GAMMA);
                interface_filename = sprintf('%s/interfaces/interface_%d.txt', ...
                    parameter_dir, k);
                displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
                
                
                
                % Loads in interface points
                [start_points, end_points] = read_interface_points(interface_filename, false); 
                
                % Finds unique values of y in all the points
                all_points = [start_points; end_points];
                [~, uniq_idxs, ~] = uniquetol(all_points(:, 2), 1e-4);
                uniq_points = all_points(uniq_idxs, :);
                
                % Sorts the resulting vector in increasing y order
                [~, sorted_idxs] = sort(uniq_points(:, 2));
                sorted_points = uniq_points(sorted_idxs, :);
                
                % Plots interface
                plot(sorted_points(:, 1), sorted_points(:, 2), 'linewidth', 2, 'displayname', displayname);
                hold on;
                
                % Find turnover point to put on graph
                dns_mat = dlmread(sprintf("%s/turnover_points.txt", parameter_dir));
                d_idx = find(dns_mat(:, 1) == k);
                d_x = dns_mat(d_idx, 2);
                d_y = dns_mat(d_idx, 3);
%                 scatter(d_x, d_y);
                
                

                
            end
        end
    end
    grid on;
    xlim([0.5, 0.55]);
    ylim([0, 0.1]);
    xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
    ylabel("$z$", "interpreter", "latex", "Fontsize", 18);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
    legend("interpreter", "latex", "location", "northwest");
    title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
    
    
    %% Figure settings
    x0=400;
    y0=400;
    width=1200;
    height=400;

    set(gcf,'position',[x0,y0,width,height])
    drawnow;
    
    
    
    
    
end
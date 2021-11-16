clear;
close all;

addpath("normal_modes");
addpath("interface_analysis");

%% Physical parameters
EPSILON = 1;
ALPHA = 1;
BETA = 10;
GAMMA = 2;
L = 16;

%% Dimensional variables
V = 5;
R = 1e-3;

millimetre_scale = R / 1e-3;
millisecond_scale = (R / V) / 1e-3;

%% Time variables
IMPACT_TIME = 0.125;
T_MAX = 0.4;
DELTA_T = 1e-4;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
TS_ANALYTICAL = 0 : DELTA_T : T_MAX - IMPACT_TIME;
IMPACT_TIMESTEP = IMPACT_TIME / DELTA_T;
NO_TIMESTEPS = length(T_VALS);

%% Finite difference parameters
N_MEMBRANE = 10924;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';


%% Data definitions
analytical_master_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping/beta_varying";
analytical_parameter_dir = sprintf("%s/alpha_%g-beta_%d-gamma_%g", analytical_master_dir, ALPHA, BETA, GAMMA);
composite_dir = sprintf("%s/finite_differences/composite", analytical_parameter_dir);
outer_dir = sprintf("%s/finite_differences/outer", analytical_parameter_dir);
nm_dir = sprintf("%s/normal_modes", analytical_parameter_dir);

dns_master_dir = "/media/michael/newarre/elastic_membrane/basilisk_parameter_sweeping/tension_varying";
dns_parameter_dir = sprintf("%s/alpha_%g-beta_%d-gamma_%g/raw_data", dns_master_dir, ALPHA, BETA, GAMMA);

stationary_dns_dir = "/media/michael/newarre/elastic_membrane/stationary_membrane/raw_data";



%% Load in DNS turnover point data
turnover_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", dns_parameter_dir));
ds_dns = turnover_mat(:, 2);

turnover_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", stationary_dns_dir));
ds_stationary_dns = turnover_mat(:, 2);

%% Load in normal modes solution
N_mat = matfile(sprintf("%s/N.mat", nm_dir));
N_nm = N_mat.N;

as_mat = matfile(sprintf("%s/as.mat", nm_dir));
as = as_mat.as;

a_ts_mat = matfile(sprintf("%s/a_ts.mat", nm_dir));
a_ts = a_ts_mat.a_ts;

q_ts_mat = matfile(sprintf("%s/q_ts.mat", nm_dir));
q_ts = q_ts_mat.q_ts;

nm_mat = matfile(sprintf("%s/ds.mat", nm_dir));
ds_nm = nm_mat.ds;

%% Figure setting up
fontsize = 20;

% Figure 1 setting up
close(figure(1));
fig1 = figure(1);
layout = tiledlayout(3, 1);
layout.Padding = 'compact';

nexttile(1, [2, 1]);
stationary_line = animatedline('color', [255, 73, 73] / 255, 'linewidth', 4);
moving_line = animatedline('color', [0 0 0], 'linewidth', 4);
upper_d_line = animatedline('color', [0 0 0], 'linewidth', 2, 'linestyle', '--');

xlim([0, 1.25]);
ylim([0, 0.5]);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', fontsize);
xticks([0 : 0.1 : 1.2]);
set(gca,'Xticklabel',[]) 
% xlabel("$x^*$ (mm)", "Fontsize", fontsize, "Interpreter", "latex");
ylabel("$z^*$ (mm)", "Fontsize", fontsize, "Interpreter", "latex");
legend(["DNS: Rigid membrane", "DNS: Deformable membrane", "Turnover point"], ...
    'Location', 'Northwest', "Interpreter", "latex")

% Figure 2 setting up
nexttile(3, [1, 1]);
ws_dns_line = animatedline('color', [0 0 0], 'linewidth', 4);
ws_fd_line = animatedline('color', 0.5 * [1 1 1], 'linestyle', '--', 'linewidth', 4);
ws_nm_line = animatedline('color', 0.5 * [1 1 1], 'linestyle', ':', 'linewidth', 4);
lower_d_line = animatedline('color', [0 0 0], 'linewidth', 2, 'linestyle', '--');

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', fontsize);
xlim([0, 1.25]);
xlabel("$x^*$ (mm)", "Fontsize", fontsize, "Interpreter", "latex");
xticks([0 : 0.1 : 1.2]);
ylabel("$-w^*(x^*, t^*)$ (mm)", "Fontsize", fontsize, "Interpreter", "latex");
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat(gca, '%.5f');
legend([ws_dns_line, ws_fd_line, ws_nm_line], ["DNS", "Analytical: Finite differences", "Analytical: Normal modes"], ...
    'Location', 'southeast', "Interpreter", "latex")


x0=400;
y0=100;
width=1800;
height=801;
set(gcf,'position',[x0,y0,width,height])

%% Video writers set up
interface_obj = VideoWriter('interfaces.avi');
open(interface_obj);
membrane_obj = VideoWriter('membrane.avi');
open(membrane_obj);

%% Loop over time
for k = 1100 : 10 : length(T_VALS)
    % Updates time
    t = T_VALS(k)
    
    %% DNS interface plotting
    nexttile(1, [2, 1]);
    title(sprintf("$t^*$ = %.4f (ms)", t * millisecond_scale), "Interpreter", "latex");
    
    % Adds interfacial lines
    stationary_filename = sprintf("%s/interface_%d.txt", stationary_dns_dir, k);
    moving_filename = sprintf("%s/interface_%d.txt", dns_parameter_dir, k);
    for filename = [moving_filename, stationary_filename, ]
        [start_points, end_points] = read_interface_points(filename, true); 
        all_points = [start_points; end_points];
        [~, uniq_idxs, ~] = uniquetol(all_points(:, 1), 1e-4);
        uniq_points = all_points(uniq_idxs, :);
        [~, sorted_idxs] = sort(uniq_points(:, 1));
        sorted_points = uniq_points(sorted_idxs, :);
        x_interp = @(z) interp1(sorted_points(:, 1), sorted_points(:, 2), z);
        
        
        if (k <= 1250)
            x_vals = zeros(1e4, 1);
            z_vals = zeros(1e4, 1);
            z_vals(2 : end) = linspace(min(sorted_points(:, 1)), max(sorted_points(:, 1)), 1e4 - 1);
            x_vals(2 : end) = x_interp(z_vals(2 : end));
        
            z_vals(1) = min(z_vals(2:end));
        else
            z_vals = linspace(min(sorted_points(:, 1)), max(sorted_points(:, 1)), 1e4);
            x_vals = x_interp(z_vals);
        end
        
        
        
        if (filename == stationary_filename)
            clearpoints(stationary_line);
            addpoints(stationary_line, x_vals, z_vals);
        else
            clearpoints(moving_line);
            addpoints(moving_line, x_vals, z_vals);
        end
    end
    hold off;
    
    if (ds_dns(k) > 0)
        clearpoints(upper_d_line);
        addpoints(upper_d_line, ds_dns(k) * [1 1], [0 1]);
    end
    
    
    
    %% Membrane plotting
    nexttile(3, [1, 1]);
    dns_ws_mat = readmatrix(sprintf("%s/w_%d.txt", dns_parameter_dir, k - 1));
    unsorted_xs_dns = dns_ws_mat(:, 1);
    [xs_dns, sort_idxs] = sort(unsorted_xs_dns);
    ws_dns = dns_ws_mat(sort_idxs, 2);
    
    clearpoints(ws_dns_line);
    addpoints(ws_dns_line, xs_dns, -ws_dns);
    
    if (t >= 0)
        %% Analytical membrane plotting
        
        % Composite plot
        ws_mat = matfile(sprintf("%s/w_%d.mat", composite_dir, k - IMPACT_TIMESTEP));
        if (k - IMPACT_TIMESTEP == 1)
            ws_composite = EPSILON^2 * ws_mat.w;
        else
            ws_composite = EPSILON^2 * ws_mat.w_next;
        end
        clearpoints(ws_fd_line);
        addpoints(ws_fd_line, xs, -ws_composite);
        
%         % Outer plot
%         ws_mat = matfile(sprintf("%s/w_%d.mat", outer_dir, k - IMPACT_TIMESTEP));
%         if (k - IMPACT_TIMESTEP == 1)
%             ws_outer = EPSILON^2 * ws_mat.w;
%         else
%             ws_outer= EPSILON^2 * ws_mat.w_next;
%         end
%         plot(xs, -ws_outer);
        
        % Normal modes plot
        [ws_nm, w_ts_nm, ps_nm] ...
            = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
            a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
            ds_nm(k - IMPACT_TIMESTEP), L, N_nm, EPSILON); 
        
        clearpoints(ws_nm_line);
        addpoints(ws_nm_line, xs, -ws_nm);
        
    end
    if (ds_dns(k) > 0)
        clearpoints(lower_d_line);
        addpoints(lower_d_line, ds_dns(k) * [1 1], [-1 1]);
    end
    
    if (t <= 0)
        ylim([min(-1.1 * max(ws_dns), -1.5e-4), max(0, -1.1 * min(ws_dns))]);
    else
        ylim([-1.1 * max([max(ws_nm), max(ws_composite), max(ws_dns)]), max([0, -1.5 * min(ws_nm)])]);
    end

    
    drawnow;
    
    f1 = getframe(fig1);
    writeVideo(interface_obj, f1);
    
    pause(0.01);
end
close(interface_obj);


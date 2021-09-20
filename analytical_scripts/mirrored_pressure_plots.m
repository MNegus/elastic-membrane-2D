%% plot_solutions.m
% Plots the saved solutions using normal modes, FD and DNS
clear;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

% Options (set to 0 if don't want to plot the solution)
normal_modes = 1;
finite_differences_comp = 1;
finite_differences_outer = 0;
dns = 1;


%% Parameters
EPSILON = 1;
ALPHA = 1 / EPSILON^2; 
BETA = 0.4 * EPSILON^2; 
GAMMA = 0.4 * EPSILON^2; 
L = 16;
T_MAX = 0.4;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/confirmation_data/beta_varying/";
analytical_parent_dir = sprintf("%s/analytical_data/alpha_%g_beta_%g_gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
dns_dir = sprintf("%s/basilisk_data/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
stationary_dns_dir = "/media/michael/newarre/elastic_membrane/confirmation_data/stationary_benchmark";

%% Loads in normal modes solutions
if (normal_modes)
    N_mat = matfile(sprintf("%s/normal_modes/N.mat", analytical_parent_dir));
    N = N_mat.N
    
    as_mat = matfile(sprintf("%s/normal_modes/as.mat", analytical_parent_dir));
    as = as_mat.as;

    a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", analytical_parent_dir));
    a_ts = a_ts_mat.a_ts;

    q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", analytical_parent_dir));
    q_ts = q_ts_mat.q_ts;
    
    % Stationary q_ts
    stationary_q_ts = stationary_normal_modes_solution(ts_analytical, N, L, EPSILON);
end



%% Turnover point compare
% Stationary dns
stationary_dns_mat = dlmread(sprintf("%s/turnover_points.txt", stationary_dns_dir));
stationary_ds_dns = stationary_dns_mat(:, 2);
stationary_ts_dns = (stationary_dns_mat(:, 1) - IMPACT_TIMESTEP) * DELTA_T;

% Find a polynomial representation of ds for positive time
idx = find(stationary_dns_mat(:, 1) == IMPACT_TIMESTEP);
stationary_ds_pos = stationary_ds_dns(idx : end);
stationary_ds_dns_sq = stationary_ds_pos.^2;
poly_stationary = polyfit(ts_analytical, stationary_ds_dns_sq, 100);
ds_dns_poly_stationary = sqrt(polyval(poly_stationary, ts_analytical)); 

% Differentiate the polynomial to get d'(t)
poly_der_stationary = polyder(poly_stationary);
stationary_d_ts_dns = (1 ./ (2 *  stationary_ds_pos)) .* polyval(poly_der_stationary, ts_analytical)';
%
if (dns)
    % DNS turnover points
    dns_mat = dlmread(sprintf("%s/turnover_points.txt", dns_dir));
    ds_dns = dns_mat(:, 2);
    ts_dns = (dns_mat(:, 1) - IMPACT_TIMESTEP) * DELTA_T;
    
    % Find a polynomial representation of ds for positive time
    idx = find(dns_mat(:, 1) == IMPACT_TIMESTEP);
    ds_pos = ds_dns(idx : end);
    ds_dns_sq = ds_pos.^2;
    poly = polyfit(ts_analytical, ds_dns_sq, 100);
    ds_dns_poly = sqrt(polyval(poly, ts_analytical)); 
    
    % Differentiate the polynomial to get d'(t)
    p_der = polyder(poly);
    d_ts_dns = (1 ./ (2 *  ds_pos)) .* polyval(p_der, ts_analytical)';
    
end

if (normal_modes)
    % Normal modes turnover points
    nm_mat = matfile(sprintf("%s/normal_modes/ds.mat", analytical_parent_dir));
    ds_nm = nm_mat.ds;
end

if (finite_differences_comp)
    % FD turnover points
    fd_comp_mat = matfile(sprintf("%s/finite_differences/composite/ds.mat", analytical_parent_dir));
    ds_comp = fd_comp_mat.ds;
end

if (finite_differences_outer)
    fd_outer_mat = matfile(sprintf("%s/finite_differences/outer/ds.mat", analytical_parent_dir));
    ds_outer = fd_outer_mat.ds;
end

%% Turnover point derivative compare
if (normal_modes)
    % Normal modes turnover points
    d_ts_nm = zeros(size(ts_analytical));
    d_ts_nm(2 : end) = diff(ds_nm) / DELTA_T;
end

if (finite_differences_comp)
    % FD turnover points
    fd_comp_mat = matfile(sprintf("%s/finite_differences/composite/d_ts.mat", analytical_parent_dir));
    d_ts_comp = fd_comp_mat.d_ts;
end

if (finite_differences_outer)
    fd_outer_mat = matfile(sprintf("%s/finite_differences/outer/d_ts.mat", analytical_parent_dir));
    d_ts_outer = fd_outer_mat.d_ts;
end

%% Loops and plots
% Line colors
colors = [[0, 0.4470, 0.7410]; ...
    [0.8500, 0.3250, 0.0980]; ...
    [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; ...
    [0.4660, 0.6740, 0.1880]; ...
    [0.3010, 0.7450, 0.9330]; ...
    [0.6350, 0.0780, 0.1840]];

% Animated lines
close all;


% for k = IMPACT_TIMESTEP - 20 : 10 : length(T_VALS)
for k = IMPACT_TIMESTEP + [51, 501, 1001, 2001] 
    %% Updates time
    t = T_VALS(k);
    t

    %% Loads in analytical solutions
    if (normal_modes)
        % Normal modes
        [ws_nm, w_ts_nm, ps_nm] ...
            = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
            a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
            ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON); 
        
        % Stationary normal modes
        [~, ~, stationary_ps_nm] ...
            = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
            a_ts(k - IMPACT_TIMESTEP, :), stationary_q_ts(k - IMPACT_TIMESTEP, :), ...
            ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON); 
    end

    if (finite_differences_comp)
        ps_comp_mat = matfile(sprintf("%s/finite_differences/composite/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
        ps_comp = ps_comp_mat.p;
    end

    if (finite_differences_outer)
        ps_outer_mat = matfile(sprintf("%s/finite_differences/outer/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
        ps_outer = ps_outer_mat.p;
    end
        
    
    %% Pressure plot
    subplot(3, 1, 3);
    xlim([0, 2]);
    
    % Stationary DNS
    stationary_pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", stationary_dns_dir, k - 1));
    unsorted_xs = stationary_pressure_mat(:, 1);
    unsorted_ps = stationary_pressure_mat(:, 2);
    [sorted_xs, idxs] = sort(unsorted_xs);
    stationary_ps = unsorted_ps(idxs);

    clearpoints(stationary_dns_p_line);
    addpoints(stationary_dns_p_line, sorted_xs, stationary_ps);
    
    % DNS
    if (dns)
        pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", dns_dir, k - 1));
        unsorted_xs = pressure_mat(:, 1);
        unsorted_ps = pressure_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        ps = unsorted_ps(idxs);
        
        clearpoints(dns_p_line);
        addpoints(dns_p_line, sorted_xs, ps);
    end
    
    
    if (normal_modes)
        clearpoints(nm_p_line);
        addpoints(nm_p_line, xs, ps_nm);
    end

    if (finite_differences_comp)
        clearpoints(comp_p_line);
        addpoints(comp_p_line, xs, ps_comp);

    end

    if (finite_differences_outer)
        clearpoints(outer_p_line);
        addpoints(outer_p_line, xs, ps_outer);
    end

    % Stationary pressure line
    clearpoints(stationary_comp_p_line);
    addpoints(stationary_comp_p_line, xs, composite_pressure_stationary(xs, t, EPSILON));
% 
%     xlim([0, 0.1]);
% 
%     if (t > 0)
%         clearpoints(d_line_3);
%         addpoints(d_line_3, ds_comp(k - IMPACT_TIMESTEP) * ones(2, 1), [-100, 100]);
% %         xline(ds_comp(k - IMPACT_TIMESTEP), 'linestyle', '--', 'linewidth', 2);
%         ylim([0, 5 * ps_comp(1)]);
%     end

    legend( "interpreter", "latex");
    xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
    ylabel("$p(x, t)$", "interpreter", "latex", "Fontsize", 18);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);

    %% Figure settings
    x0=400;
    y0=400;
    width=1200;
    height=800;

    set(gcf,'position',[x0,y0,width,height])
    drawnow;
    frame = getframe(gcf);
%         writeVideo(writerobj, frame);

    pause(2);
end
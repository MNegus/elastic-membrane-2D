%% plot_solutions.m
% Plots the saved solutions using normal modes, FD and DNS
clear;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

% Options (set to 0 if don't want to plot the solution)
normal_modes = 1;
finite_differences_comp = 0;
finite_differences_outer = 0;
dns = 0;



%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters()
ALPHA = ALPHAS(1);
BETA = BETAS(1);
GAMMA = GAMMAS(1);

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX;

%% Data dirs
parent_dir = "/home/michael/Desktop/alpha_vary_test";
analytical_parent_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
dns_dir = "/media/michael/newarre/elastic_membrane/gamma_vary_test/basilisk_data/gamma_0.1";

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
end

%% Turnover point compare

% if (dns)
%     % DNS turnover points
%     ds_dir = "/home/michael/scratch/coupled_benchmark";
%     dns_mat = dlmread(sprintf("%s/turnover_points.txt", ds_dir));
%     ds_dns = dns_mat(:, 2);
%     ts_dns = (dns_mat(:, 1) - IMPACT_TIMESTEP) * DELTA_T;
% end

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

close(figure(1));
figure(1);
hold on;
% plot(ts_dns, ds_dns, 'linewidth', 2);
plot(ts_analytical, ds_nm, 'linewidth', 2);
% plot(ts_analytical, ds_outer, 'linewidth', 2);
% plot(ts_analytical, ds_comp, 'linewidth', 2);
legend(["DNS", "Normal modes", "FD: Outer", "FD: Composite"], 'location', 'northwest');
xlabel("t");
ylabel("d(t)");

%% Loops and plots
% writerobj = VideoWriter("four_model_compare.avi");
% writerobj.FrameRate = 10;
% open(writerobj);

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

if (normal_modes)
    subplot(3, 1, 1);
    nm_w_line = animatedline('Color', colors(3, :),'LineWidth', 2, ...
        'Displayname', 'Normal modes');
    
    subplot(3, 1, 2);
    nm_w_t_line = animatedline('Color', colors(3, :),'LineWidth', 2, ...
        'Displayname', 'Normal modes');
    
    subplot(3, 1, 3);
    nm_p_line = animatedline('Color', colors(3, :),'LineWidth', 2, ...
        'Displayname', 'Normal modes');
end

if (finite_differences_comp)
    subplot(3, 1, 1);
    comp_w_line = animatedline('Color', colors(2, :),'LineWidth', 2, ...
        'Displayname', 'FD: Composite');
    
    subplot(3, 1, 2);
    comp_w_t_line = animatedline('Color', colors(2, :),'LineWidth', 2, ...
        'Displayname', 'FD: Composite');
    
    subplot(3, 1, 3);
    comp_p_line = animatedline('Color', colors(2, :),'LineWidth', 2, ...
        'Displayname', 'FD: Composite');
end

if (finite_differences_outer)
    subplot(3, 1, 1);
    outer_w_line = animatedline('Color', colors(4, :),'LineWidth', 2, ...
        'Displayname', 'FD: Outer');
    
    subplot(3, 1, 2);
    outer_w_t_line = animatedline('Color', colors(4, :),'LineWidth', 2, ...
        'Displayname', 'FD: Outer');
    
    subplot(3, 1, 3);
    outer_p_line = animatedline('Color', colors(4, :),'LineWidth', 2, ...
        'Displayname', 'FD: Outer');
    
end


if (dns)
    subplot(3, 1, 1);
    dns_w_line = animatedline('Color', colors(1, :),'LineWidth', 2, ...
        'Displayname', 'DNS');
    
    subplot(3, 1, 2);
    dns_w_t_line = animatedline('Color', colors(1, :),'LineWidth', 2, ...
        'Displayname', 'DNS');
    
    subplot(3, 1, 3)
    dns_p_line = animatedline('Color', colors(1, :),'LineWidth', 2, ...
        'Displayname', 'DNS');
end




% d line
subplot(3, 1, 1);
d_line_1 = animatedline('LineWidth', 2, 'linestyle', '--', ...
        'Displayname', 'Turnover point');
subplot(3, 1, 2);
d_line_2 = animatedline('LineWidth', 2, 'linestyle', '--', ...
        'Displayname', 'Turnover point');
subplot(3, 1, 3);
d_line_3 = animatedline('LineWidth', 2, 'linestyle', '--', ...
        'Displayname', 'Turnover point');


for k = IMPACT_TIMESTEP - 200 : 100 : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t

    %% Loads in analytical solutions
    if (t > 0)
        
        if (normal_modes)
            % Normal modes
            [ws_nm, w_ts_nm, ps_nm] ...
                = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
                a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), ...
                ds_nm(k - IMPACT_TIMESTEP), L, N, EPSILON); 
        end
        
        if (finite_differences_comp)
            % Composites
            ws_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ws_comp = EPSILON^2 * ws_comp_mat.w_next;

            w_ts_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            w_ts_comp = EPSILON^2 * w_ts_comp_mat.w_t;

            ps_comp_mat = matfile(sprintf("%s/finite_differences/composite/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ps_comp = ps_comp_mat.p;
        end
        
        if (finite_differences_outer)
            % Outers
            ws_outer_mat = matfile(sprintf("%s/finite_differences/outer/w_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ws_outer = ws_outer_mat.w_next;

            w_ts_outer_mat = matfile(sprintf("%s/finite_differences/outer/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            w_ts_outer = w_ts_outer_mat.w_t;

            ps_outer_mat = matfile(sprintf("%s/finite_differences/outer/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ps_outer = ps_outer_mat.p;
        end
    end
        
    %% w plot
    subplot(3, 1, 1);
%     set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');

    % DNS
    if (dns)
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", dns_dir, k - 1));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        
        clearpoints(dns_w_line);
        addpoints(dns_w_line, sorted_xs, ws);
    end

    if (t > 0)
        if (normal_modes)
            clearpoints(nm_w_line);
            addpoints(nm_w_line, xs, ws_nm);
        end
        
        if (finite_differences_comp)
            clearpoints(comp_w_line);
            addpoints(comp_w_line, xs, ws_comp);
        end
        
        if (finite_differences_outer)
            clearpoints(outer_w_line);
            addpoints(outer_w_line, xs, ws_outer);
        end
    end
    
%     if (t > 0)
%             clearpoints(d_line_1);
%             addpoints(d_line_1, ds_comp(k - IMPACT_TIMESTEP) * ones(2, 1), [-100, 100]);
%     %         xline(ds_comp(k - IMPACT_TIMESTEP), 'linestyle', '--', 'linewidth', 2);
%             ylim([0, 1.2 * ws_comp(1)]);
%         end
    % 

    xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
    ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
    legend("interpreter", "latex");


    title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

    %% w_t plot
    subplot(3, 1, 2);
    
    % DNS
    if (dns)
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", dns_dir, k - 1));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_w_ts = membrane_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        w_ts = unsorted_w_ts(idxs);
        
        clearpoints(dns_w_t_line);
        addpoints(dns_w_t_line, sorted_xs, w_ts);
        
    end

    if (t > 0)
        if (normal_modes)
            clearpoints(nm_w_t_line);
            addpoints(nm_w_t_line, xs, w_ts_nm);
        end
        
        if (finite_differences_comp)
            clearpoints(comp_w_t_line);
            addpoints(comp_w_t_line, xs, w_ts_comp);
        end
        
        if (finite_differences_outer)
            clearpoints(outer_w_t_line);
            addpoints(outer_w_t_line, xs, w_ts_outer);
        end
    end
    
%     if (t > 0)
%         clearpoints(d_line_2);
%         addpoints(d_line_2, ds_comp(k - IMPACT_TIMESTEP) * ones(2, 1), [-100, 100]);
% %         xline(ds_comp(k - IMPACT_TIMESTEP), 'linestyle', '--', 'linewidth', 2);
%         ylim([0, 1.2 * w_ts_comp(1)]);
%     end
% 
    xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
    ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", 18);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
    legend("interpreter", "latex");

    %% Pressure plot
    subplot(3, 1, 3);
    
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
    

    if (t > 0)
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
    end
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

    pause(0.1);
end
% close(writerobj);
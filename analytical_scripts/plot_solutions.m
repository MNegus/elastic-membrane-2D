%% save_solutions.m
% Plots the saved solutions using normal modes, FD and DNS
clear;

addpath("finite_differences");
addpath("normal_modes");
addpath("pressures");

analytical_parent_dir = "/media/michael/newarre/elastic_membrane/analytical_tests";
dns_dir = "~/scratch/reflecting_waves/membrane_radius_4";

%% Parameters
EPSILON = 1;
ALPHA = 2 / EPSILON^2; BETA = 1 * EPSILON^2; GAMMA = 2 * EPSILON^2; 
L = 4;
T_MAX = 0.25;
DELTA_T = 1e-4;


% NM parameters
N = 128;

% FD parameters
N_MEMBRANE = 2056;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX;

%% Loads in normal modes solutions
as_mat = matfile(sprintf("%s/normal_modes/as.mat", analytical_parent_dir));
as = as_mat.as;

a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", analytical_parent_dir));
a_ts = a_ts_mat.a_ts;

q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", analytical_parent_dir));
q_ts = q_ts_mat.q_ts;

%% Turnover point compare
ts_analytical = 0 : DELTA_T : T_MAX;

% DNS turnover points
ds_dir = "/home/michael/scratch/coupled_benchmark";
dns_mat = dlmread(sprintf("%s/turnover_points.txt", ds_dir));
ds_dns = dns_mat(:, 2);
ts_dns = (dns_mat(:, 1) - IMPACT_TIMESTEP) * DELTA_T;

% Normal modes turnover points
nm_mat = matfile(sprintf("%s/normal_modes/ds.mat", analytical_parent_dir));
ds_nm = nm_mat.ds;

% FD turnover points
fd_comp_mat = matfile(sprintf("%s/finite_differences/composite/ds.mat", analytical_parent_dir));
ds_comp = fd_comp_mat.ds_composite;

fd_outer_mat = matfile(sprintf("%s/finite_differences/outer/ds.mat", analytical_parent_dir));
ds_outer = fd_outer_mat.ds_outer;

close(figure(1));
figure(1);
hold on;
plot(ts_dns, ds_dns, 'linewidth', 2);
plot(ts_analytical, ds_nm, 'linewidth', 2);
plot(ts_analytical, ds_outer, 'linewidth', 2);
plot(ts_analytical, ds_comp, 'linewidth', 2);
legend(["DNS", "Normal modes", "FD: Outer", "FD: Composite"], 'location', 'northwest');
xlabel("t");
ylabel("d(t)");

%% Loops and plots
% writerobj = VideoWriter("four_model_compare.avi");
% writerobj.FrameRate = 10;
% open(writerobj);

% for k = IMPACT_TIMESTEP - 200 : length(T_VALS)
for k = 1500
    %% Updates time
    t = T_VALS(k);
    t

    %% Plots
%     if ((mod(k-1, 10) == 0))
        
        %% Loads in analytical solutions
        if (t > 0)
            
            % Normal modes
            [ws_nm, w_ts_nm, ps_nm] ...
                = w_solution_normal_modes(xs, as(k - IMPACT_TIMESTEP, :), ...
                a_ts(k - IMPACT_TIMESTEP, :), q_ts(k - IMPACT_TIMESTEP, :), L, N); 
            
            % Composites
            ws_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ws_comp = ws_comp_mat.w_next_composite;
            
            w_ts_comp_mat = matfile(sprintf("%s/finite_differences/composite/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            w_ts_comp = w_ts_comp_mat.w_t_composite;
            
            ps_comp_mat = matfile(sprintf("%s/finite_differences/composite/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ps_comp = ps_comp_mat.p_composite;

            % Outers
            ws_outer_mat = matfile(sprintf("%s/finite_differences/outer/w_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ws_outer = ws_outer_mat.w_next_outer;
            
            w_ts_outer_mat = matfile(sprintf("%s/finite_differences/outer/w_t_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            w_ts_outer = w_ts_outer_mat.w_t_outer;
            
            ps_outer_mat = matfile(sprintf("%s/finite_differences/outer/p_%d.mat", analytical_parent_dir, k - IMPACT_TIMESTEP));
            ps_outer = ps_outer_mat.p_outer;
        end
        
        %% w plot
        subplot(3, 1, 1);
        set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
        
        % Reads in Basilisk solution
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", dns_dir, k - 1));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        plot(sorted_xs, ws, 'linewidth', 2);
        
        if (t > 0)
            hold on;
            plot(xs, ws_nm, 'linewidth', 2);
            plot(xs, ws_outer, 'linewidth', 2);
            plot(xs, ws_comp, 'linewidth', 2);
            hold off;
        end
        
        xlim([0, 0.1]);
        
        if (t > 0)
            xline(ds_comp(k - IMPACT_TIMESTEP), 'linestyle', '--', 'linewidth', 2);
            xlim([0, max(0.1, 5 * ds_comp(k - IMPACT_TIMESTEP))]);
            
        end

        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        legend(["DNS", "Normal modes", "FD: Outer", "FD: Composite"], "interpreter", "latex");


        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

        %% w_t plot
        subplot(3, 1, 2);
        % Reads in Basilisk solution
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", dns_dir, k - 1));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_w_ts = membrane_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        w_ts = unsorted_w_ts(idxs);
        plot(sorted_xs, w_ts, 'linewidth', 2);
        
        if (t > 0)
            hold on;
            plot(xs, w_ts_nm, 'linewidth', 2);
            plot(xs, w_ts_outer, 'linewidth', 2);
            plot(xs, w_ts_comp, 'linewidth', 2);
            hold off;
        end

        xlim([0, 0.1]);
        
        if (t > 0)
            xline(ds_comp(k - IMPACT_TIMESTEP), 'linestyle', '--', 'linewidth', 2);
            xlim([0, max(0.1, 5 * ds_comp(k - IMPACT_TIMESTEP))]);
        end
        
        legend(["DNS", "Normal modes", "FD: Outer", "FD: Composite", "Turnover point"], "interpreter", "latex");
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        legend(["DNS", "Normal modes", "FD: Outer", "FD: Composite"], "interpreter", "latex");

        %% Pressure plot
        subplot(3, 1, 3);
%         figure(7);
        % Reads in Basilisk solution
        pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", dns_dir, k - 1));
        unsorted_xs = pressure_mat(:, 1);
        unsorted_ps = pressure_mat(:, 2);
        [sorted_xs, idxs] = sort(unsorted_xs);
        ps = unsorted_ps(idxs);
        plot(sorted_xs, ps, 'linewidth', 2);
        
        if (t > 0)
            hold on;
            plot(xs, ps_nm, 'linewidth', 2);
            plot(xs, ps_outer, 'linewidth', 2);
            plot(xs, ps_comp, 'linewidth', 2);
            hold off;
        end

        xlim([0, 0.1]);
        
        if (t > 0)
            xline(ds_comp(k - IMPACT_TIMESTEP), 'linestyle', '--', 'linewidth', 2);
            xlim([0, max(0.1, 5 * ds_comp(k - IMPACT_TIMESTEP))]);
            ylim([0, 5 * ps_comp(1)]);
        end
        
        legend(["DNS", "Normal modes", "FD: Outer", "FD: Composite", "Turnover point"], "interpreter", "latex");
        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$p(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);

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
% end

% close(writerobj);
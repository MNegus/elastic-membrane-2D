%% save_solutions.m
% Plots the saved solutions using normal modes
clear;

analytical_parent_dir = "~/Desktop/analytical_tests";

%% Parameters
EPSILON = 1;
ALPHA = 80 / EPSILON^2; BETA = 0 * EPSILON^2; GAMMA = 1e6 * EPSILON^2; 
L = 4;
T_MAX = 0.25;
% T_MAX = 5e-4;
DELTA_T = 1e-4;

% NM parameters
N = 32;

N_MEMBRANE = 1024;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';


%% Loads in normal modes solutions
as_mat = matfile(sprintf("%s/normal_modes/as.mat", analytical_parent_dir));
as = as_mat.as;

a_ts_mat = matfile(sprintf("%s/normal_modes/a_ts.mat", analytical_parent_dir));
a_ts = a_ts_mat.a_ts;

q_ts_mat = matfile(sprintf("%s/normal_modes/q_ts.mat", analytical_parent_dir));
q_ts = q_ts_mat.q_ts;

%% Turnover point compare
ts_analytical = 0 : DELTA_T : T_MAX;

% Normal modes turnover points
nm_mat = matfile(sprintf("%s/normal_modes/ds.mat", analytical_parent_dir));
ds_nm = nm_mat.ds;

close(figure(2));
figure(2);
hold on;
plot(ts_analytical, 2 * sqrt(ts_analytical), 'linewidth', 4);
plot(ts_analytical, ds_nm, 'linewidth', 2);
legend(["Stationary", "Normal modes"], 'location', 'northwest');
xlabel("t");
ylabel("d(t)");

%% Loops and plots
% writerobj = VideoWriter("four_model_compare.avi");
% writerobj.FrameRate = 10;
% open(writerobj);

for k = 1 : 100 : length(ts_analytical)
    %% Updates time
    t = ts_analytical(k);
    t

    %% Plots
%     if ((mod(k-1, 10) == 0))
        
        %% Loads in analytical solutions
        % Normal modes
        [ws_nm, w_ts_nm, ps_nm] ...
            = w_solution_normal_modes(xs, as(k, :), ...
            a_ts(k, :), q_ts(k, :), L, N); 
        
        %% w plot
        subplot(3, 1, 1);
        set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
        
        plot(xs, ws_nm, 'linewidth', 2);
        xline(ds_nm(k), 'linestyle', '--', 'linewidth', 2);
%         xlim([0, max(0.1, 5 * ds_nm(k))]);

        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        legend(["Normal modes"], "interpreter", "latex");


        title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 

        %% w_t plot
        subplot(3, 1, 2);
        set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
        
        plot(xs, w_ts_nm, 'linewidth', 2);
        xline(ds_nm(k), 'linestyle', '--', 'linewidth', 2);
%         xlim([0, max(0.1, 5 * ds_nm(k))]);

        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        legend(["Normal modes"], "interpreter", "latex");

        %% Pressure plot
        subplot(3, 1, 3);
       set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
        
        plot(xs, ps_nm, 'linewidth', 2);
        xline(ds_nm(k), 'linestyle', '--', 'linewidth', 2);
%         xlim([0, max(0.1, 5 * ds_nm(k))]);

        xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
        ylabel("$p(x, t)$", "interpreter", "latex", "Fontsize", 18);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
        legend(["Normal modes"], "interpreter", "latex");

        %% Figure outputting
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
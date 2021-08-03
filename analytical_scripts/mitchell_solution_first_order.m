%% mitchell_solution_stationary.m
% Code to use an iterative Mitchell FD method

clear;
close all;

addpath("pressures");
addpath("finite_differences");

%% Parameters
EPSILON = 1;
ALPHA = 2 / EPSILON^2; BETA = 1 * EPSILON^2; GAMMA = 2 * EPSILON^2; 
L = 16;
N_MEMBRANE = 10924;
M = N_MEMBRANE - 1; % We ignore the end point
T_MAX = 0.25;
DELTA_T = 10^-4;
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = IMPACT_TIME / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX;

% Derived parameters
if (GAMMA == 0)
    Cpressure = DELTA_X * DELTA_X / BETA;
else
    Cpressure = DELTA_X^4 / GAMMA;
end

%% Basilisk compare data
data_directory = "~/scratch/reflecting_waves/membrane_radius_16";

%% Matrix definitions
[A_mat, B_mat] = matrix_definitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);

%% Initialise arrays

% Composite solutions
w_previous_composite = zeros(size(xs));
w_composite = initialise_membrane(w_previous_composite, A_mat, B_mat);
w_t_composite = zeros(size(xs));
w_next_composite = zeros(size(xs));
p_previous_previous_composite = zeros(size(xs));
p_previous_composite = zeros(size(xs));
p_composite = zeros(size(xs));

% Outer solutions
w_previous_outer = zeros(size(xs));
w_outer = initialise_membrane(w_previous_outer, A_mat, B_mat);
w_t_outer = zeros(size(xs));
w_next_outer = zeros(size(xs));
p_previous_previous_outer = zeros(size(xs));
p_previous_outer = zeros(size(xs));
p_outer = zeros(size(xs));

%% Stationary solutions to time-dependent terms
d_fun = @(t) 2 * sqrt(t);
d_t_fun = @(t) 1 / sqrt(t);
A_fun = @(t) d_fun(t) * d_t_fun(t);
C_fun = @(t) d_fun(t) * d_t_fun(t);
J_fun = @(t) pi * d_fun(t) / (8 * d_t_fun(t)^2);
    
%% Loops
writerobj = VideoWriter("first_order_numerical_comparison.avi");
writerobj.FrameRate = 10;
open(writerobj);

for k = IMPACT_TIMESTEP : length(T_VALS)
    %% Updates time
    t = T_VALS(k);
    t
    
    %% Analytical solution
    if k > IMPACT_TIMESTEP
        
        %% Composite timestep
        tic
        [w_next_composite, p_composite, w_t_composite, d, d_t] = membrane_timestep(xs, t, ...
            w_composite, w_previous_composite, p_previous_composite, p_previous_previous_composite, ...
            "composite",  ALPHA, BETA, GAMMA, EPSILON, ...
            M, DELTA_X, DELTA_T, Cpressure, A_mat, B_mat);

        % Swaps ws
        temp = w_previous_composite;
        w_previous_composite = w_composite;
        w_composite = w_next_composite;
        w_next_composite = temp;

        % Swaps ps
        temp = p_previous_previous_composite;
        p_previous_previous_composite = p_previous_composite;
        p_previous_composite = p_composite;
        p_composite = temp;
        toc
        %% Outer timestep
        [w_next_outer, p_outer, w_t_outer, d, d_t] = membrane_timestep(xs, t, ...
            w_outer, w_previous_outer, p_previous_outer, p_previous_previous_outer, ...
            "outer",  ALPHA, BETA, GAMMA, EPSILON, ...
            M, DELTA_X, DELTA_T, Cpressure, A_mat, B_mat);

        % Swaps ws
        temp = w_previous_outer;
        w_previous_outer = w_outer;
        w_outer = w_next_outer;
        w_next_outer = temp;

        % Swaps ps
        temp = p_previous_previous_outer;
        p_previous_previous_outer = p_previous_outer;
        p_previous_outer = p_outer;
        p_outer = temp;
        
    end

%     %% Plots
%     if (mod(k-1, 10) == 0)
%         % w plot
%         subplot(3, 1, 1);
%         set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
%         % Reads in Basilisk solution
%         membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", data_directory, k - 1));
%         unsorted_xs = membrane_mat(:, 1);
%         unsorted_ws = membrane_mat(:, 2);
%         [sorted_xs, idxs] = sort(unsorted_xs);
%         ws = unsorted_ws(idxs);
%         plot(sorted_xs, ws, 'linewidth', 2);
% 
%         if (t > 0)
%             hold on;
%             plot(xs, w_next_composite, 'linewidth', 2);
%             plot(xs, w_next_outer, 'linewidth', 2);
%             hold off;
%         end
% 
%         xlim([0, 4]);
%         xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
%         ylabel("$w(x, t)$", "interpreter", "latex", "Fontsize", 18);
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%         legend(["DNS", "Analytical: Composite", "Analytical: Outer"], "interpreter", "latex");
% 
% 
%         title(sprintf("$t$ = %.4f", t), "Interpreter", "latex"); 
% 
%         % w_t plot
%         subplot(3, 1, 2);
%         % Reads in Basilisk solution
%         membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", data_directory, k - 1));
%         unsorted_xs = membrane_mat(:, 1);
%         unsorted_w_ts = membrane_mat(:, 2);
%         [sorted_xs, idxs] = sort(unsorted_xs);
%         w_ts = unsorted_w_ts(idxs);
%         plot(sorted_xs, w_ts, 'linewidth', 2);
%         if (t > 0)
%             hold on;
%             plot(xs, w_t_composite, 'linewidth', 2);
%             plot(xs, w_t_outer, 'linewidth', 2);
%             hold off;
%         end
%         xlim([0, 4]);
%         xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
%         ylabel("$w_t(x, t)$", "interpreter", "latex", "Fontsize", 18);
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%         legend(["DNS", "Analytical: Composite", "Analytical: Outer"], "interpreter", "latex");
% 
%         % Pressure plot
%         subplot(3, 1, 3);
%         % Reads in Basilisk solution
%         pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", data_directory, k - 1));
%         unsorted_xs = pressure_mat(:, 1);
%         unsorted_ps = pressure_mat(:, 2);
%         [sorted_xs, idxs] = sort(unsorted_xs);
%         ps = unsorted_ps(idxs);
%         plot(sorted_xs, ps, 'linewidth', 2);
% 
%         if (t > 0)
%             hold on;
%             plot(xs, p_previous_composite, 'linewidth', 2);
%             plot(xs, p_previous_outer, 'linewidth', 2);
% 
%             p_stationary = composite_pressure_stationary(xs, t, d_fun(t), d_t_fun(t), A_fun(t), C_fun(t), J_fun(t), EPSILON);
%             plot(xs, p_stationary, 'linewidth', 2);
%             
%             hold off;
%             ylim([0, 5 * p_composite(1)]);
%         end
%         legend(["DNS", "Analytical: Composite", "Analytical: Outer", "Analytical: Stationary"], "interpreter", "latex");
% 
%         xlim([0, 4]);
%         xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
%         ylabel("$p(x, t)$", "interpreter", "latex", "Fontsize", 18);
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%         
%         
% 
%         x0=400;
%         y0=400;
%         width=1200;
%         height=800;
% 
%         set(gcf,'position',[x0,y0,width,height])
% 
%         frame = getframe(gcf);
%         writeVideo(writerobj, frame);
%         pause(0.00001);
%     end
end
close(writerobj);

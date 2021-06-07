%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution

clear;
close all;
%% Parameters

L = 4;
N_MEMBRANES = [1365, 2730, 5461, 10922];
N_legend = strings(size(N_MEMBRANES));
for n = 1 : length(N_legend)
   N_legend(n) = "$N$ = " + N_MEMBRANES(n);
end

Deltaxs =  L ./ (N_MEMBRANES - 1)

T_MAX = 0.4;
DELTA_POWERS = linspace(-1, -5, 9)
DELTA_TS = 10.^DELTA_POWERS
DT_legend = strings(size(DELTA_TS));
for k = 1 : length(DT_legend)
   DT_legend(k) = "$\Delta t$ = " + DELTA_TS(k); 
end
DT_legend


ALPHA = 1;
BETA = 1;
GAMMA = 1;
% 
% Dbetas = BETA * DELTA_TS.^2 ./ (ALPHA * Deltaxs'.^2)
% Dgammas = GAMMA * DELTA_TS.^2 ./ (ALPHA * Deltaxs'.^4)



%% Numerically solves and compare
max_errors = zeros(length(DELTA_TS), length(Deltaxs));
L2_norms = zeros(length(DELTA_TS), length(Deltaxs));


for n = 1 : length(N_MEMBRANES)
    N_MEMBRANE = N_MEMBRANES(n);
    M = N_MEMBRANE - 1; % We ignore the end point
    xs = linspace(0, L, N_MEMBRANE);
    Deltax = L / (N_MEMBRANE - 1)

    for k = 1 : length(DELTA_TS)
        % Reset figure
%         close(figure(1));
%         figure(1);
%         
%         close(figure(2));
%         figure(2);

        % Setting DELTA_T
        DELTA_T = DELTA_TS(k);
        Dbeta = BETA * DELTA_T^2 / (ALPHA * Deltax^2)
        Dgamma = GAMMA * DELTA_T^2 / (ALPHA * Deltax^4)

        % Initialises error term
        max_error = 0;

        % A definition
        A_upper_upper = Dgamma * ones(M, 1);
        A_upper_upper(3) = 2 * Dgamma;

        A_upper = (-Dbeta - 4 * Dgamma) * ones(M, 1);
        A_upper(2) = -2 * Dbeta - 8 * Dgamma;

        A_main = (4 + 2 * Dbeta + 6 * Dgamma) * ones(M, 1);
        A_main(2) = A_main(2) + Dgamma;
        A_main(M) = A_main(M) - Dgamma;

        A_lower = (-Dbeta - 4 * Dgamma) * ones(M, 1);

        A_lower_lower = Dgamma * ones(M, 1);

        A = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);

        % B definition
        B_upper_upper = -2 * Dgamma * ones(M, 1);
        B_upper_upper(3) = 2 * B_upper_upper(3);

        B_upper = (2 * Dbeta + 8 * Dgamma) * ones(M, 1);
        B_upper(2) = 2 * B_upper(2);

        B_main = (8 - 4 * Dbeta - 12 * Dgamma) * ones(M, 1);
        B_main(2) = B_main(2) - 2 * Dgamma;
        B_main(M) = B_main(M) + 2 * Dgamma;

        B_lower = (2 * Dbeta + 8 * Dgamma) * ones(M, 1);

        B_lower_lower = -2 * Dgamma * ones(M, 1);

        B = spdiags([B_lower_lower B_lower B_main B_upper B_upper_upper], -2:2, M, M);

        % Initial conditions
        t = 0
        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
        w_previous = w_exact(1 : end - 1);

%         plot(xs(1 : end - 1), w_previous);
%         hold on;
%         plot(xs, w_exact);
%         hold off;
%         pause(0.01);

        % Solve for first w
        t = t + DELTA_T
        rhs = 0.5 * B * w_previous;
        w = A \ rhs;

        w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);


        % Initialise error norm
        L2_norm = sum((w - w_exact(1 : end - 1)).^2);

%         Loops
        max_diff = 0
        while (t < T_MAX) 
            rhs = B * w - A * w_previous;
            w_next = A \ rhs;

            % Swaps
            temp = w_previous;
            w_previous = w;
            w = w_next;
            w_next = temp;

            t = t + DELTA_T

            % Exact solution
            w_exact = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);

            % Updates max errors 
            max_diff = max(abs(w - w_exact(1 : end - 1)));
            if (max_diff > max_error)
               max_error = max_diff; 
            end
            
            L2_norm = L2_norm + sum((w - w_exact(1 : end - 1)).^2);
            
%             % Plots both
%             figure(1);
%             plot(xs(1 : end - 1), w);
%             hold on;
%             plot(xs, w_exact);
%             hold off;
%             
%             % Plots diff
%             figure(2);
%             plot(xs(1 : end - 1), w - w_exact(1 : end - 1));
%             pause(0.01);

        end
        max_errors(k, n) = max_diff;
        
        % Saves L2_norm
        L2_norm = sqrt(L2_norm / (length(xs) * (T_MAX / DELTA_T)));
        L2_norms(k, n) = L2_norm
    end
    
%     figure(3);
%     hold on;
%     plot(DELTA_TS, max_errors, '-o');
%     
%     figure(4);
%     hold on;
%     plot(DELTA_TS, L2_norms, '-o');
    
end

%% Max norm vs DT
close(figure(3));
figure(3);
hold on;
for n = 1 : length(N_MEMBRANES)
    plot(DELTA_TS, max_errors(:, n), '-o', 'linewidth', 1.5);
end
plot(DELTA_TS, 0.3 * DELTA_TS.^2, 'linestyle', '--', 'color', 'black', 'linewidth', 2);
set(gca, 'XDir','reverse');
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize', 12);
% set(gca,'ticklabel
legend([N_legend, "$y \sim \Delta t^2$"], "Interpreter", "latex", "fontsize", 12);
xlabel("$\Delta t$", "interpreter", "latex", "fontsize", 12);
ylabel("Max norm error", "interpreter", "latex", "fontsize", 12);
title(["$\alpha$ = " + ALPHA + ", $\beta$ = " + BETA + ", $\gamma$ = " + GAMMA + ", $L$ = " + L], "Interpreter", "latex");
exportgraphics(gca, "timestep_validation_max.png", "resolution", 300);

%% Max norm vs DX
close(figure(4));
figure(4);
hold on;
for k = 1 : length(DELTA_TS)
    plot(N_MEMBRANES, max_errors(k, :), '-o', 'linewidth', 1.5);
end
plot(N_MEMBRANES, 0.1./ N_MEMBRANES.^2, 'linestyle', '--', 'color', 'black', 'linewidth', 2);
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize', 12);
legend([DT_legend, "$y \sim 1 / N^2$"], "Interpreter", "latex", "fontsize", 12, "location", "eastoutside");
xlabel("$N$", "interpreter", "latex", "fontsize", 12);
ylabel("Max norm error", "interpreter", "latex", "fontsize", 12);
title(["$\alpha$ = " + ALPHA + ", $\beta$ = " + BETA + ", $\gamma$ = " + GAMMA + ", $L$ = " + L], "Interpreter", "latex");
exportgraphics(gca, "gridsize_validation_max.png", "resolution", 300);


%% L2 norm vs DT
close(figure(5));
figure(5);
hold on;
for n = 1 : length(N_MEMBRANES)
    plot(DELTA_TS, L2_norms(:, n), '-o', 'linewidth', 1.5);
end
plot(DELTA_TS, 0.1 * DELTA_TS.^2, 'linestyle', '--', 'color', 'black', 'linewidth', 2);
set(gca, 'XDir','reverse');
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize', 12);
% set(gca,'ticklabel
legend([N_legend, "$y \sim \Delta t^2$"], "Interpreter", "latex", "fontsize", 12);
xlabel("$\Delta t$", "interpreter", "latex", "fontsize", 12);
ylabel("$L_2$ norm error", "interpreter", "latex", "fontsize", 12);
title(["$\alpha$ = " + ALPHA + ", $\beta$ = " + BETA + ", $\gamma$ = " + GAMMA + ", $L$ = " + L], "Interpreter", "latex");
exportgraphics(gca, "timestep_validation_L2.png", "resolution", 300);

%% L2 norm vs DX
close(figure(6));
figure(6);
hold on;
for k = 1 : length(DELTA_TS)
    plot(N_MEMBRANES, L2_norms(k, :), '-o', 'linewidth', 1.5);
end
plot(N_MEMBRANES, 0.03./ N_MEMBRANES.^2, 'linestyle', '--', 'color', 'black', 'linewidth', 2);
set(gca, 'yscale','log');
set(gca, 'xscale', 'log');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize', 12);
legend([DT_legend, "$y \sim 1 / N^2$"], "Interpreter", "latex", "fontsize", 12, "location", "eastoutside");
xlabel("$N$", "interpreter", "latex", "fontsize", 12);
ylabel("$L_2$ norm error", "interpreter", "latex", "fontsize", 12);
title(["$\alpha$ = " + ALPHA + ", $\beta$ = " + BETA + ", $\gamma$ = " + GAMMA + ", $L$ = " + L], "Interpreter", "latex");
exportgraphics(gca, "gridsize_validation_L2.png", "resolution", 300);

% %% L2 norm plot 
% figure(4);
% set(gca, 'XDir','reverse');
% set(gca, 'yscale','log');
% set(gca, 'xscale', 'log');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize', 12);
% xlabel("$\Delta t$", "interpreter", "latex", "fontsize", 12);
% ylabel("$L_2$ norm error", "interpreter", "latex", "fontsize", 12);
% legend(["$N$ = 512", "$N$ = 1024", "$N$ = 2048", "$N$ = 4096", "$N$ = 8192"], "Interpreter", "latex", "fontsize", 12);
% title(["$\alpha$ = " + ALPHA + ", $\beta$ = " + BETA + ", $\gamma$ = " + GAMMA + ", $L$ = " + L], "Interpreter", "latex");
% exportgraphics(gca, "timestep_validation_L2.png");



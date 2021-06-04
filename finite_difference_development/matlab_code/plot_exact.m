%% mitchell_method_comparison.m
% Code to compare the Mitchell method scheme for the wave equation solution
% to the exact solution
close all;
clear;
%% Parameters
L = 4;
N_MEMBRANE = 8192;
Deltaxs =  L ./ (N_MEMBRANE - 1)
T_MAX = 0.4;
DELTA_T = 1e-3;

% ALPHA = 1;
% BETA = ALPHA * (L / (5461 - 1))^2 / ((1e-4)^2)
% GAMMA = ALPHA * (L / (5461 - 1))^4 / ((1e-4)^2)
ALPHA = 1; BETA = 1; GAMMA = 1;


tvals = 0 : DELTA_T : T_MAX;
xs = linspace(0, L, N_MEMBRANE);

%% Numerically solves and compare
writerobj = VideoWriter("exact_video.avi");
open(writerobj);

figure(1);
for t = tvals
    t
    % Exact solution
    w_exact_1 = homogeneous_exact_solution(xs, t, ALPHA, BETA, GAMMA, L);
%     w_exact_2 = homogeneous_exact_solution(xs, t, ALPHA, BETA, 0, L);

    %  Plots full
    plot(xs, w_exact_1, 'linewidth', 2);
    ylim([-0.02, 0.18]);
    xlim([0, L]);
    fontsize = 20;
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'fontsize', fontsize);
    xlabel("$x$", "interpreter", "latex", "fontsize", fontsize);
    ylabel("$w(x, t)$", "interpreter", "latex", "fontsize", fontsize);
    grid on;

    frame = getframe(gcf);
    writeVideo(writerobj,frame);
    pause(0.01);

end
close(writerobj);



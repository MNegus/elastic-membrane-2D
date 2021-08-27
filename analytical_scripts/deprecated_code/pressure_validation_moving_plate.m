%% Pressure validation, moving plate
% Validates the computations for the pressure, assuming that w(x, t) = w(t)
% where w(t) = 0.5 + t - 0.5 * sqrt(1 + 4 * t).
% close all;
addpath("~/chebfun");
addpath("pressures");

%% Parameters
L = 1;
epsilon = 0.1;
N = 1024;
dx = L / (N - 1);
xs = (0 : dx : L)';
tmax = (L / (2 * epsilon))^2;
dt = tmax / 100;
tvals = dt : dt : tmax;


%% w definitions (moving plate)
w_fun = @(x, t) (0.5 + t - 0.5 * sqrt(1 + 4 * t)) * ones(size(x));
w_t_fun = @(x, t) (1 - 1 / sqrt(1 + 4 * t)) * ones(size(x));
w_tt_fun = @(x, t) (1 / (2 * (1 + 4 * t)^1.5)) * ones(size(x));
w_x_fun = @(x, t) zeros(size(x));
m_t_fun = @(s, t) w_t_fun(epsilon * s, t) .* s;
m_tt_fun = @(s, t) w_tt_fun(epsilon * s, t) .* s;

%% w definitions (zeros)
% w_fun = @(x, t) zeros(size(x));
% w_t_fun = @(x, t) zeros(size(x));
% w_tt_fun = @(x, t) zeros(size(x));
% w_x_fun = @(x, t) zeros(size(x));
% m_t_fun = @(x, t) zeros(size(x));
% m_tt_fun = @(x, t) zeros(size(x));


%% Arrays for saved quantities
exact_ds = zeros(size(tvals));
exact_d_ts = zeros(size(tvals));
exact_As = zeros(size(tvals));
exact_Cs = zeros(size(tvals));
exact_Js = zeros(size(tvals));

numerical_ds = zeros(size(tvals));
numerical_d_ts = zeros(size(tvals));
numerical_As = zeros(size(tvals));
numerical_Cs = zeros(size(tvals));
numerical_Js = zeros(size(tvals));

%% Timings arrays
time_dependents_times = zeros(size(tvals));
outer_times = zeros(size(tvals));
inner_times = zeros(size(tvals));
overlap_times = zeros(size(tvals));

%% Loop over time
writerobj = VideoWriter("moving_plate_pressure.avi");
writerobj.FrameRate = 10;
open(writerobj);
figure(1);
for k = 1 : length(tvals)
    t = tvals(k)

    % Anonymous functions for w quantities
    w = @(x) w_fun(x, t);
    w_t = @(x) w_t_fun(x, t);
    w_tt = @(x) w_tt_fun(x, t);
    w_x = @(x) w_x_fun(x, t);
    m_t = @(x) m_t_fun(x, t);
    m_tt = @(x) m_tt_fun(x, t);

    % Determine numerical solutions
    tic;
    [d, d_t] = turnover_point(t, w, w_t, w_x, epsilon)
    [A, C, J] = time_dependent_quantities(d, d_t, w_t, w_tt, m_t, epsilon)
    time_dependents_times(k) = toc;
    
    % Save numerical solutions
    numerical_ds(k) = d;
    numerical_d_ts(k) = d_t;
    numerical_As(k) = A;
    numerical_Cs(k) = C;
    numerical_Js(k) = J;
   
    % Determine exact solutions
    exact_ds(k) = 2 * sqrt(t - w(0));
    exact_d_ts(k) = (1 - w_t(0)) / sqrt(t - w(0));
    exact_As(k) = d * d_t * (1 - w_t(0)) - w_tt(0) * d^2 / 2;
    exact_Cs(k) = d * d_t * (1 - w_t(0));
    exact_Js(k) = pi * (1 - w_t(0))^2 * d / (8 * d_t^2);

    % Determine pressure solutions
    tic;
    outer_ps = outer_pressure(xs, m_tt, w_tt, d, A, epsilon);
    outer_times(k) = toc;
    
    tic;
    overlap_ps = overlap_pressure(xs, d, C, epsilon);
    overlap_times(k) = toc;
    
    tic;
    inner_ps = inner_pressure(xs, d, d_t, J, epsilon);
    inner_times(k) = toc;
    
    composite_ps = outer_ps + inner_ps - overlap_ps;
    
    pmax = max(composite_ps);
   
    figure(1);
    plot(xs, outer_ps, 'linestyle', '--', 'linewidth', 2, 'color', 0.75 * [1 1 1]);
    hold on;
    plot(xs, inner_ps, 'linestyle', ':', 'linewidth', 2, 'color', 0.75 * [1 1 1]);
    plot(xs, composite_ps, 'linewidth', 2, 'color', 'black');
    hold off;
    legend("Outer", "Inner", "Composite", "Interpreter", "latex", "fontsize", 15, ...
        "location", "northeast");
    ylim([-0.001, 1.2 * pmax]);
    xlabel("$x$", "Interpreter", "latex", "fontsize", 20);
	ylabel("$p(x, t)$", "Interpreter", "latex", "fontsize", 20);
    title("$t = $" + t, "Interpreter", "latex", "fontsize", 20);
    
    x0=400;
    y0=400;
    width=1000;
    height=500;
    
    set(gcf,'position',[x0,y0,width,height])
    frame = getframe(gcf);
    writeVideo(writerobj, frame);
end
close(writerobj);

%% Plot times compare
figure(2);
plot(tvals, time_dependents_times, 'linewidth', 2);
hold on;
plot(tvals, outer_times, 'linewidth', 2);
plot(tvals, inner_times, 'linewidth', 2);
plot(tvals, overlap_times, 'linewidth', 2);
hold off;
set(gca, "yscale", "log");
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("Computational time / s", "Interpreter", "latex", "fontsize", 20);
legend(["$d(t)$, $d'(t)$, etc calculation", "Outer calculation", "Inner calculation", "Overlap calculation"], "Interpreter", "latex", "fontsize", 12, ...
        "location", "east");
exportgraphics(gcf, "overall_times.png");

%% Plot d and d_t solutions
subplot(2, 1, 1);
d_norms = abs((exact_ds - numerical_ds) ./ exact_ds);
plot(tvals, d_norms);
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("$||d_0(t)||$", "Interpreter", "latex", "fontsize", 20);
xlim([-1e-3, tmax]);
ylim([-0.1 * max(d_norms), 1.2 * max(d_norms)]);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 15);
grid on;

subplot(2, 1, 2);
d_t_norms = abs((exact_d_ts - numerical_d_ts) ./ exact_d_ts);
plot(tvals, d_t_norms);
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("$||d'_0(t)||$", "Interpreter", "latex", "fontsize", 20);
xlim([-1e-3, tmax]);
ylim([-0.1 * max(d_t_norms), 1.2 * max(d_t_norms)]);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 15);
grid on;

x0=400;
y0=400;
width=1000;
height=500;
set(gcf,'position',[x0,y0,width,height]);
sgtitle("Turnover point validation: Moving plate case", "Interpreter", "latex", "fontsize", 20) 
exportgraphics(gcf, "turnover_validation_moving_plate.png");

%% A, C, J validation
subplot(3, 1, 1);
A_norms = abs((exact_As - numerical_As) ./ exact_As);
plot(tvals, A_norms);
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("$||A(t)||$", "Interpreter", "latex", "fontsize", 20);
xlim([-1e-3, tmax]);
ylim([-0.1 * max(A_norms), 1.2 * max(A_norms)]);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 15);
grid on;



subplot(3, 1, 2);
C_norms = abs((exact_Cs - numerical_Cs) ./ exact_Cs);
plot(tvals, C_norms);
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("$||C(t)||$", "Interpreter", "latex", "fontsize", 20);
xlim([-1e-3, tmax]);
ylim([-0.1 * max(C_norms), 1.2 * max(C_norms)]);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 15);
grid on;

subplot(3, 1, 3);
J_norms = abs((exact_Js - numerical_Js) ./ exact_Js);
plot(tvals, J_norms);
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("$||J(t)||$", "Interpreter", "latex", "fontsize", 20);
xlim([-1e-3, tmax]);
ylim([-0.1 * max(J_norms), 1.2 * max(J_norms)]);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 15);
grid on;

x0=400;
y0=400;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height]);
sgtitle("A, B, C, J validation: Moving plate case", "Interpreter", "latex", "fontsize", 20) 
exportgraphics(gcf, "ABCJ_validation_moving_plate.png");

%% Pressure validation, moving plate
% Validates the computations for the pressure, assuming that w(x, t) = w(t)
% where w(t) = 0.5 + t - 0.5 * sqrt(1 + 4 * t).
% close all;
addpath("pressures");

%% Parameters
L = 1;
epsilon = 0.1;
N = 1024;
dx = L / (N - 1);
xs = 0 : dx : L;
tmax = (L / (2 * epsilon))^2;
dt = tmax / 100;
tvals = dt : dt : tmax;

% L = 1;
% epsilon = 1;
% N = 1024;
% dx = L / (N - 1);
% xs = 0 : dx : L;
% tmax = (L / (2 * epsilon))^2;
% dt = tmax / 1000;
% tvals = dt : dt : tmax;

%% w definitions (moving plate)
w_fun = @(x, t) (0.5 + t - 0.5 * sqrt(1 + 4 * t)) * ones(size(x));
w_t_fun = @(x, t) (1 - 1 / sqrt(1 + 4 * t)) * ones(size(x));
w_tt_fun = @(x, t) (1 / (2 * (1 + 4 * t)^1.5)) * ones(size(x));
w_x_fun = @(x, t) zeros(size(x));
m_t_fun = @(s, t) w_t_fun(epsilon * s, t) .* s;
m_tt_fun = @(s, t) w_tt_fun(epsilon * s, t) .* s;


%% Arrays for saved quantities
exact_ds = zeros(size(tvals));
exact_d_ts = zeros(size(tvals));
exact_As = zeros(size(tvals));
exact_Bs = zeros(size(tvals));
exact_Cs = zeros(size(tvals));
exact_Js = zeros(size(tvals));

numerical_ds = zeros(size(tvals));
numerical_d_ts = zeros(size(tvals));
numerical_As = zeros(size(tvals));
numerical_Bs = zeros(size(tvals));
numerical_Cs = zeros(size(tvals));
numerical_Js = zeros(size(tvals));

%% Loop over time
figure(1);
integral_times = zeros(size(tvals));
trapz_times = zeros(size(tvals));

writerobj = VideoWriter("outer_pressure_compare.avi");
writerobj.FrameRate = 10;
open(writerobj);
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
    
    [d, d_t, A, B, C, J] ...
       = time_dependent_quantities(t, w, w_t, w_tt, w_x, m_t, epsilon);
    
    % Determine pressure solutions
    tic
    integral_outer_ps = outer_pressure_integral_method(xs, m_tt, d, A, epsilon);
    toc
    integral_times(k) = toc;
    
    
    tic
    trapz_outer_ps = outer_pressure(xs, m_tt, w_tt, d, A, epsilon);
    toc
    trapz_times(k) = toc;
    
    
    figure(1);
    plot(xs, integral_outer_ps, 'linewidth', 5, 'linestyle', '--', 'color', 0.5 * [1 1 1]);
    hold on;
    plot(xs, trapz_outer_ps, 'linewidth', 2, 'color', 'black');
    hold off;
    legend(["integral outer", "trapz outer"], "Interpreter", "latex", "fontsize", 15, ...
        "location", "northeast");
    xlabel("$x$", "Interpreter", "latex", "fontsize", 20);
	ylabel("$p(x, t)$", "Interpreter", "latex", "fontsize", 20);
    
    tstring = sprintf("%.2f", t);
    title("$t = $" + tstring, "Interpreter", "latex", "fontsize", 20);
    ylim([0, integral_outer_ps(1) * 10]);
    
    x0=400;
    y0=400;
    width=1000;
    height=500;
    
    set(gcf,'position',[x0,y0,width,height]);
    
    set(gcf,'position',[x0,y0,width,height])
    frame = getframe(gcf);
    writeVideo(writerobj, frame);
end
close(writerobj);

%% times plot
figure(2);
plot(tvals, integral_times);
hold on;
plot(tvals, trapz_times);
hold off;
set(gca, "yscale", "log");
xlabel("$t$", "Interpreter", "latex", "fontsize", 20);
ylabel("Computational time / s", "Interpreter", "latex", "fontsize", 20);
legend(["integral method", "trapz method"], "Interpreter", "latex", "fontsize", 15, ...
        "location", "east");
exportgraphics(gcf, "integral_trapz_times.png");
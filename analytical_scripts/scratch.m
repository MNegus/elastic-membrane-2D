addpath("pressures");

clear;
close all;

%% Parameters
L = 1;
epsilon = 1;
N = 1024;
dx = L / (N - 1);
xs = 0 : dx : L;

%% w definition

% Stationary
% w_fun = @(s) zeros(size(s));
% w_t_fun = @(s) zeros(size(s));
% w_tt_fun = @(s) zeros(size(s));
% m_tt_fun = @(s) zeros(size(s));

% Constant velocity
% w_fun = @(s) t^2 * ones(size(s));
% w_t_fun = @(s) 2 * t * ones(size(s));
% w_tt_fun = @(s) 2 * ones(size(s));
% m_tt_fun = @(s) 2 * s;

% Oscillating mode
% omega = 2 * pi;
% a = 0.1
% w_fun = @(s) a * sin(omega * t) * cos(lambda * s);
% w_x_fun = @(s) - a * lambda * sin(omega * t) * sin(lambda * s);
% w_t_fun = @(s) a * omega * cos(omega * t) * cos(lambda * s);
% w_tt_fun = @(s) - a * omega^2 * sin(omega * t) * cos(lambda * s);
% m_t_fun = @(s) a * omega * cos(omega * t) * sin(lambda * s) / lambda;
% m_tt_fun = @(s) - a * omega^2 * sin(omega * t) * sin(lambda * s) / lambda;

%%

for t = 1e-9 : 1e-3 : 0.25
% for t = 0.1
    
    %% w definitions 
    % Oscillating mode
    omega = 2 * pi;
    lambda = pi / (2 * L);
    a = 0.1;
    w_fun = @(s) a * sin(omega * t) * cos(lambda * s);
    w_x_fun = @(s) - a * lambda * sin(omega * t) * sin(lambda * s);
    w_t_fun = @(s) a * omega * cos(omega * t) * cos(lambda * s);
    w_tt_fun = @(s) - a * omega^2 * sin(omega * t) * cos(lambda * s);
    m_t_fun = @(s) a * omega * cos(omega * t) * sin(lambda * s) / lambda;
    m_tt_fun = @(s) - a * omega^2 * sin(omega * t) * sin(lambda * s) / lambda;
    
    %% Determine time-dependent quantities
    
    % Determine d(t)
    d_zero_fun = @(d) full_d_zero_fun(d, t, w_fun, epsilon);
    d = fsolve(d_zero_fun, 2 * sqrt(t));
    
    % Determine d'(t)
    d_t_zero_fun = @(d_t) full_d_t_zero_fun(d_t, d, w_t_fun, w_x_fun, epsilon);
    d_t = fsolve(d_t_zero_fun, 1 / sqrt(t));

    % Determine J(t)
    B_integrand = @(s) (m_t_fun(d) - m_t_fun(s)) ./ (sqrt(d^2 - s.^2) .* (d - s));
    B = (1 / pi) * integral(B_integrand, -d, d);
    J = pi * (1 - B)^2 * d / (8 * d_t^2);

    outer_ps = outer_pressure(xs, w_t_fun, w_tt_fun, m_tt_fun, d, d_t, epsilon);
    overlap_ps = overlap_pressure(xs, w_t_fun, d, d_t, epsilon);
    inner_ps = inner_pressure(xs, d, d_t, J, epsilon);
    composite_ps = outer_ps + inner_ps - overlap_ps;
    figure(1);
    plot(xs, outer_ps, 'linestyle', '--', 'linewidth', 2);
    hold on;
    plot(xs, inner_ps, 'linestyle', ':', 'linewidth', 2);
    plot(xs, composite_ps, 'linewidth', 2);
    hold off;
    legend("Outer", "Inner", "Composite");
%     hold on;
%     plot(xs, 2 ./(epsilon * sqrt(epsilon^2 * d^2 - xs.^2)));
%     hold off;
%     ylim([0, 1.25 * max(composite_ps)]);
    ylim([0, 50]);
    title("t = " + t);
    pause(0.01);
end

%% Compare to exact
% integral_1_exact = 2 * t * besselj(0, epsilon * lambda * d);
% integral_1_exact - integral_1
% 
% integral_2_exact = 2 * d * besselj(1, epsilon * lambda * d) / (epsilon * lambda);
% integral_2_exact - integral_2


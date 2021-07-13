addpath("pressures");

close all;

L = 1;
epsilon = 0.1;
N = 256;
dx = L / (N - 1);
xs = 0 : dx : L;

% E.g. w(x, t) = t^2 cos(lambda x)
lambda = pi / (2 * L);

for t = 0 : 1e-3 : 1
    t
    c = 10 * 2 * pi;
%     w_fun = @(s) sin(c * t) * cos(lambda * s);
%     w_t_fun = @(s) c * cos(t) * cos(lambda * s);
%     w_tt_fun = @(s) - c^2 * sin(t) * cos(lambda * s);
%     m_tt_fun = @(s) - c^2 * sin(t) * sin(lambda * s) / lambda;

    w_fun = @(s) zeros(size(s));
    w_t_fun = @(s) zeros(size(s));
    w_tt_fun = @(s) zeros(size(s));
    m_tt_fun = @(s) zeros(size(s));

    % w_fun = @(s) t^2 * ones(size(s));
    % w_t_fun = @(s) 2 * t * ones(size(s));
    % w_tt_fun = @(s) 2 * ones(size(s));
    % m_tt_fun = @(s) 2 * s;


    d = 2 * sqrt(t);
    d_t = 1 / sqrt(t);

    outer_ps = outer_pressure(xs, w_t_fun, w_tt_fun, m_tt_fun, d, d_t, epsilon);
    overlap_ps = overlap_pressure(xs, w_t_fun, d, d_t, epsilon);
    figure(1);
    plot(xs, outer_ps);
    hold on;
    plot(xs, -overlap_ps);
    plot(xs, outer_ps - overlap_ps);
    hold off;
%     hold on;
%     plot(xs, 2 ./(epsilon * sqrt(epsilon^2 * d^2 - xs.^2)));
%     hold off;
%     ylim([0, 500]);
    pause(0.01);
end

%% Compare to exact
% integral_1_exact = 2 * t * besselj(0, epsilon * lambda * d);
% integral_1_exact - integral_1
% 
% integral_2_exact = 2 * d * besselj(1, epsilon * lambda * d) / (epsilon * lambda);
% integral_2_exact - integral_2
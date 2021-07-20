xs = linspace(0, 1, 1e4);
epsilon = 0.1;


% zero_eqn = @(eta) (2 * eta + 1 - 4 * exp(-eta) - exp(-2 * eta))
% eta_0 = fsolve(zero_eqn, 0);



% [etas, xs_etas] = etas_definition(xs, d, J, epsilon);
% 
% figure(1);
% plot(xs_etas, etas, '-o');

tvals = linspace(0.1, 10, 1e3);
d_t_vals = 1 ./ sqrt(tvals);
p_max = zeros(size(tvals));

for k = 1 : length(tvals)
    t = tvals(k)
    d = 2 * sqrt(t);
    d_t = 1 / sqrt(t);
    J = pi * d / (8 * d_t^2);
    
    ps = inner_pressure(xs, d, d_t, J, epsilon);
    p_max(k) = max(ps);
%     ps_alt = inner_pressure_alt(xs, d, d_t, J, epsilon);
    
    figure(2);
    plot(xs, ps);
%     hold on;
%     plot(xs_etas, ps_etas);
%     hold off;
    ylim([0, 50]);
end

%% 
close(figure(3));
figure(3);
plot(tvals, 2 * d_t_vals.^2 / (epsilon^2 * 4));
hold on;
plot(tvals, p_max);

legend("Theory", "P_max_xs", "P_max_etas");
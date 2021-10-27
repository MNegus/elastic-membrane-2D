clear;

%% Analytical functions
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);

%% Computational values
no_points = 512;
x_vels = zeros(no_points, 1);
y_vels = zeros(no_points, 1);
options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);
timesteps = 1260 : 10 : 3000;
impact_time = 0.125;
DELTA_T = 1e-4;
tvals = DELTA_T * timesteps - impact_time;

%% Basilisk data 
data_dir = "/scratch/negus/no_surface_condition/raw_data";
turnover_point_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", data_dir));
basilisk_d_ts = turnover_point_mat(:, 4);

%% Arrays for fluxes
fluxes_times = [];
analytical_fluxes = [];
basilisk_fluxes = [];

%% Determine etas
etas = jet_root_etas(no_points);

%% Plots velocity at jet-root at each time
close(figure(1));
figure(1);
hold on;

close(figure(2));
figure(2);
hold on;

close(figure(3));
figure(3);
analytical_flux_line = animatedline('color', 'b');
basilisk_flux_line = animatedline('color', 'r');
jet_flux_line = animatedline('color', 'g');
yline(pi / 2);
legend("Analytical", "Basilisk");

for idx = 1 : length(timesteps)
    k = timesteps(idx);
    t = tvals(idx);
    
    %% Loops over etas to find z, and the corresponding velocities
    zs = zeros(size(etas));
    
    for m = 1 : length(etas)
        eta = etas(m);
        
        % Determine z
        xz = zeta_fun(eta, J(t));
        zs(m) = xz(2);
        x = xz(1)

        % Find velocities (in vertical moving frame
        w = (1 + 1i * sqrt(eta)) / (1 - 1i * sqrt(eta));
        x_vels(m) = d_t(t) * real(w) + d_t(t);
        y_vels(m) = -d_t(t) * imag(w);
        
    end
    
    %% Reads in basilisk data
    energy_mat = readmatrix(sprintf("%s/energy_%d.txt", data_dir, k));
    basilisk_ys = energy_mat(:, 1);
    basilisk_x_vels = energy_mat(:, 3);
    basilisk_y_vels = energy_mat(:, 4);
    
    %% Determines energy fluxes
%     analytical_integrand = ((x_vels).^2 + y_vels.^2) .* (x_vels - d_t(t));
    analytical_integrand = ((x_vels).^2) .* (x_vels - d_t(t));
    
    analytical_flux = trapz(zs, analytical_integrand)
    addpoints(analytical_flux_line, t, analytical_flux);
    
    % Jet flux
    eta_jet = 0;
    w_jet = (1 + 1i * sqrt(eta_jet)) / (1 - 1i * sqrt(eta_jet));
    x_vel_jet = d_t(t) * real(w_jet) + d_t(t);
    z_jet = z_fun(eta_jet, J(t));
    jet_flux = z_jet * (x_vel_jet^2) * (x_vel_jet - d_t(t));
    addpoints(jet_flux_line, t, jet_flux);
    
    basilisk_d_t = basilisk_d_ts(k);
    basilisk_integrand = (basilisk_x_vels.^2 + basilisk_y_vels.^2) .* (basilisk_x_vels - basilisk_d_t);
    if (length(basilisk_integrand) >= 2)
        basilisk_flux = trapz(basilisk_ys, basilisk_integrand);
    else
        basilisk_flux = 0;
    end
    addpoints(basilisk_flux_line, t, basilisk_flux);
    
    
    %% Plots x velocity distribution
    figure(1);
    h(1) = plot(x_vels, zs, 'displayname', "Analytical");
    hold on;
    h(2) = plot(basilisk_x_vels, basilisk_ys, 'displayname', 'Basilisk');
    h(3) = xline(d_t(t), 'linestyle', '--');
    h(4) = yline(J(t), 'linestyle', '--');
    h(5) = yline(J(t) * (1 + 4 / pi), 'linestyle', '--');
    h(6) = xline(2 * d_t(t), 'linestyle', ':');
    legend(h(1:2));
    hold off;
    ylim([0, 1.2 * J(t) * (1 + 4 / pi)]);
    ylabel("z");
    xlabel("u_x");
    
    
    %% Plots y velocity distribution
    figure(2);
    plot(y_vels, zs);
    hold on;
    plot(basilisk_y_vels, basilisk_ys);
    hold off;
    legend();
    ylim([0, 1.2 * J(t) * (1 + 4 / pi)]);
    ylabel("y");
    xlabel("u_y");
    
    t
    drawnow;
    pause(0.01);
end


%% Function definitions
function z = z_fun(eta, J)
    z = (J / pi) * imag(eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
end

function xz = zeta_fun(eta, J)
    zeta = (J / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
    xz(1) = real(zeta);
    xz(2) = imag(zeta);
end

function eta = eta_0(options)
    eta_guess = -0.01 + 1e-10;
    eta = fsolve(@(eta) (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1), eta_guess, options);
end

% function integrand = analytical_integrand(z, t, d_t, J, options)
% 
%     eta_fun = @(eta) z - z_fun(eta, J(t));
%    
%     eta0 = 1 + 1i * 1e-6;
%     eta = fsolve(eta_fun, eta0, options);
%     
%     %% Find velocities
%     w = (1 + 1i * sqrt(eta)) / (1 - 1i * sqrt(eta));
%     x_vel = d_t(t) * real(w);
% %     y_vel = -d_t(t) * imag(w) + 1;
%     
%     %% Find integrand
%     integrand = (x_vel + d_t(t))^2 * x_vel
%     
% end
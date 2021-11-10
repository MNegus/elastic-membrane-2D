clear;

%% Analytical functions
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);

%% Computational values
no_points = 512;
u_tildes = zeros(no_points, 1);
v_tildes = zeros(no_points, 1);
p_tildes = zeros(no_points, 1);
options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);
timesteps = 1290 : 10 : 3000;
impact_time = 0.125;
DELTA_T = 1e-4;
tvals = DELTA_T * timesteps - impact_time;

%% Basilisk data 
data_dir = "/scratch/negus/flux_with_pressure/raw_data";
turnover_point_mat = readmatrix(sprintf("%s/turnover_points_basilisk.txt", data_dir));
basilisk_d_ts = turnover_point_mat(:, 4);

%% Arrays for fluxes
fluxes_times = [];
analytical_fluxes = [];
basilisk_fluxes = [];

%% Determine etas
x = 0;
etas = jet_root_etas(x, no_points);

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
legend("Analytical", "Basilisk", "Jet flux");

close(figure(4));
figure(4);
hold on;

% for idx = 1 : length(timesteps)
for idx = 1 : 100
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
        u_tildes(m) = d_t(t) * real(w);
        v_tildes(m) = -d_t(t) * imag(w);
        
        
    end
    
    p_tildes = -0.5 * (u_tildes.^2 + v_tildes.^2) + 0.5 * d_t(t)^2;
    
    %% Reads in basilisk data
    x_cell = 0;
    velocities_mat = readmatrix(sprintf("%s/velocities_%d-x_cell_%d.txt", data_dir, k, x_cell));
    basilisk_zs = velocities_mat(:, 1);
    basilisk_fs = velocities_mat(:, 2);
    basilisk_x_vels = velocities_mat(:, 3);
    basilisk_z_vels = velocities_mat(:, 4);
    basilisk_ps = velocities_mat(:, 5);
    
    %% Determines energy fluxes
    analytical_integrand = u_tildes.^2;
    analytical_flux = d_t(t)^3 * J(t) + d_t(t) * trapz(zs, analytical_integrand);

%     analytical_integrand = (0.5 * ((u_tildes + d_t(t)).^2 + v_tildes.^2) + p_tildes) .* u_tildes;
%     analytical_flux = trapz(zs, analytical_integrand);
    addpoints(analytical_flux_line, t, analytical_flux);
    
    % Jet flux
%     eta_jet = 0;
%     w_jet = (1 + 1i * sqrt(eta_jet)) / (1 - 1i * sqrt(eta_jet));
%     x_vel_jet = d_t(t) * real(w_jet) + d_t(t);
%     z_jet = z_fun(eta_jet, J(t));
%     jet_flux = z_jet * (x_vel_jet^2) * (x_vel_jet - d_t(t));
    jet_flux = 2 * d_t(t)^3 * J(t);
    addpoints(jet_flux_line, t, jet_flux);
    
    % Basilisk flux
    basilisk_d_t = basilisk_d_ts(k);
%     basilisk_integrand = (basilisk_x_vels.^2 + basilisk_y_vels.^2) .* (basilisk_x_vels - basilisk_d_t);
    basilisk_integrand = (basilisk_x_vels - basilisk_d_t) .* (0.5 * (basilisk_x_vels.^2 + basilisk_z_vels.^2) + basilisk_ps);
    if (length(basilisk_integrand) >= 2)
        basilisk_flux = trapz(basilisk_zs, basilisk_integrand);
    else
        basilisk_flux = 0;
    end
    addpoints(basilisk_flux_line, t, basilisk_flux);
    
    
    %% Plots x velocity distribution
    figure(1);
    h(1) = plot(u_tildes + d_t(t), zs, 'displayname', "Analytical");
    hold on;
    h(2) = plot(basilisk_x_vels, basilisk_zs, 'displayname', 'Basilisk');
    h(3) = xline(d_t(t), 'linestyle', '--');
    h(4) = yline(J(t), 'linestyle', '--');
    h(5) = yline(J(t) * (1 + 4 / pi), 'linestyle', '--');
    h(6) = xline(2 * d_t(t), 'linestyle', ':');
    legend(h(1:2));
    hold off;
    ylim([0, 1.2 * J(t) * (1 + 4 / pi)]);
    ylabel("z");
    xlabel("Horizontal velocity, u(x, z, t)");
    
    
    %% Plots y velocity distribution
    figure(2);
    plot(v_tildes, zs, 'Displayname', 'Analytical');
    hold on;
    plot(basilisk_z_vels, basilisk_zs, 'Displayname', 'Basilisk');
    hold off;
    legend();
    ylim([0, 1.2 * J(t) * (1 + 4 / pi)]);
    ylabel("z");
    xlabel("Vertical velocity, v(x, z, t)");
    
    
    %% Plots pressure distribution
    figure(4);
    plot(p_tildes, zs, 'Displayname', 'Analytical');
    hold on;
    plot(basilisk_ps, basilisk_zs, 'Displayname', 'Basilisk');
    hold off;
    legend();
    ylabel("z");
    xlabel("Pressure, p(x, z, t)");
    
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



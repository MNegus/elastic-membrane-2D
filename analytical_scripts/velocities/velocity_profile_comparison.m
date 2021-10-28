clear;

%% Analytical functions
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);

%% Computational values
no_points = 512;
x_vels = zeros(no_points, 1);
z_vels = zeros(no_points, 1);
options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);
timesteps = 1260 : 10 : 3000;
impact_time = 0.125;
DELTA_T = 1e-4;
tvals = DELTA_T * timesteps - impact_time;


%% x values to test profile at
x_profiles = linspace(0, 1, 4);
etas_profiles = zeros(length(x_profiles), no_points);
for q = 1 : length(x_profiles)
    etas_profiles(q, :) = jet_root_etas(x_profiles(q), no_points);
end

%% Arrays for fluxes
fluxes_times = zeros(size(timesteps));
fluxes = zeros(length(x_profiles), length(timesteps));
jet_fluxes = zeros(length(timesteps));

%% Plots velocity at jet-root at each time
close(figure(1));
figure(1);
hold on;

close(figure(2));
figure(2);
hold on;

close(figure(3));
figure(3);

for idx = 1 : length(timesteps)
    k = timesteps(idx);
    t = tvals(idx);
    flux_times(idx) = t;
    
    %% Plots free surface on each figure
    eta_heights = linspace(0, 1, 1e3);
    x_heights = (J(t) / pi) * (eta_heights - log(eta_heights) - 1);
    z_heights = J(t) * (1 + 4 * sqrt(eta_heights) / pi);
    
%     figure(1);
%     yyaxis right;
%     plot(x_heights, z_heights, 'color', 'black', 'linewidth', 2);
%     yyaxis left;
%     
%     figure(2);
%     yyaxis right;
%     plot(x_heights, z_heights, 'color', 'black', 'linewidth', 2);
%     yyaxis left;
    
    %% Loop over x_profile values, plotting their respective profiles
    for q = 1 : length(x_profiles)
        etas = etas_profiles(q, :);
        
        %% Finds xs and zs, and their velocities
        xs = zeros(size(etas)); zs = zeros(size(etas));
        x_vels = zeros(size(etas)); z_vels = zeros(size(etas));
        for m = 1 : length(etas)
            eta = etas(m);
            
            % Determine x and z
            xz = zeta_fun(eta, J(t));
            xs(m) = xz(1);
            zs(m) = xz(2);
            
            % Find velocities in inner frame
            w = (1 + 1i * sqrt(eta)) / (1 - 1i * sqrt(eta));
            x_vels(m) = d_t(t) * real(w);
            z_vels(m) = -d_t(t) * imag(w);
        end
        
        %% Determine energy fluxes 
        
        % NOT INCLUDING Y VELOCITY
        analytical_integrand = ((x_vels + d_t(t)).^2) .* (x_vels - d_t(t));
        analytical_flux = trapz(zs, analytical_integrand);
        fluxes(q, idx) = analytical_flux;
        
        %% Plot x velocity
        figure(1);
        h(q) = plot(x_vels, zs, 'displayname', sprintf("x = %g", x_profiles(q)));
        hold on;
        
        %% Plot z velocity
        figure(2);
        g(q) = plot(z_vels, zs, 'displayname', sprintf("x = %g", x_profiles(q)));
        hold on;
        
        %% Plot flux
        figure(3);
        p(q) = plot(flux_times(1 : idx), fluxes(q, 1 : idx), 'displayname', sprintf("x = %g", x_profiles(q)));
        hold on;
        
    end
    
    
    %% x velocity distribution plot
    figure(1);
    legend(h(1 : length(x_profiles)));
    ylabel("z");
    xlabel("u_x");
    hold off;
    
    %% z velocity distribution plot
    figure(2);
    legend(g(1 : length(x_profiles)));
    ylabel("z");
    xlabel("u_z");
    hold off;
    
    %% Flux distribution plot
    figure(3);
    % Jet flux
    eta_jet = 0;
    w_jet = (1 + 1i * sqrt(eta_jet)) / (1 - 1i * sqrt(eta_jet));
    x_vel_jet = d_t(t) * real(w_jet) + d_t(t);
    z_jet = z_fun(eta_jet, J(t));
    jet_fluxes(idx) = z_jet * (x_vel_jet^2) * (x_vel_jet - d_t(t));
    p(q) = plot(flux_times(1 : idx), jet_fluxes(1 : idx), 'displayname', "Jet flux");
    
    legend();
    ylabel("Flux");
    xlabel("t");
    hold off;
    
    %% Output quantities
    t
    drawnow;
    pause(1);
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

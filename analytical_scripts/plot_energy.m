%% plot_energy.m


%% Parameters
[EPSILON, ALPHAS, BETAS, GAMMAS, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX;
ts_analytical = 0 : DELTA_T : T_MAX;

% Pressure type (composite or outer)
pressure_type = "composite";


%% Data dirs
parent_dir = "/home/negus/Desktop/jet_energy_test";


%% Turnover point comparison
close(figure(1));
figure(1);
hold on;

% Plot stationary value

plot(ts_analytical, 2 * pi * ts_analytical, 'linewidth', 2, 'Displayname', "Stationary");

for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            % Loads in parameters
            parameter_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              parent_dir, ALPHA, BETA, GAMMA, pressure_type)
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];
          
            % Extract d_ts and Js
            Js_mat = matfile(sprintf("%s/Js.mat", parameter_dir));
            Js = Js_mat.Js;
            
            d_ts_mat = matfile(sprintf("%s/d_ts.mat", parameter_dir));
            d_ts = d_ts_mat.d_ts;
            
            % Determines jet energy
            jet_energy = 8 * EPSILON^3 * cumtrapz(ts_analytical, Js .* d_ts.^3);
            
            % Plot line
            plot(ts_analytical, jet_energy, 'linewidth', 2, 'Displayname', displayname);
        end
    end
end
legend("location", "northeast");
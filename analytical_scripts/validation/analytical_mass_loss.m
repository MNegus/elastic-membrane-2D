% Plots the saved solutions using normal modes, FD and DNS
clear;
% close all;

% addpath("../finite_differences");
% addpath("../normal_modes");
% addpath("../pressures");
addpath("../");
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 16;
N_MEMBRANE = 10924;
DELTA_T = 1e-4;
T_MAX = 0.4;

EPSILON = 1;

% FD parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0.125 / DELTA_T;
T_VALS = -IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

%% Stationary values
ds_stationary = 2 * sqrt(ts_analytical);
d_ts_stationary = 1 ./ sqrt(ts_analytical);
Js_stationary = pi * ds_stationary ./ (8 * d_ts_stationary.^2);
mass_flux_stationary = 2 * Js_stationary .* d_ts_stationary;
mass_flux_stationary(1) = 0;
mass_stationary = cumtrapz(ts_analytical, mass_flux_stationary);

%% Loads in moving membrane solution
moving_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/finite_differences/composite";
Js_mat = matfile(sprintf("%s/Js.mat", moving_dir));
Js = Js_mat.Js;

d_ts_mat = matfile(sprintf("%s/d_ts.mat", moving_dir));
d_ts = d_ts_mat.d_ts;

mass_flux = 2 * Js .* d_ts;
mass = cumtrapz(ts_analytical, mass_flux);
%%
figure(7)
plot(ts_analytical, mass_stationary);
hold on;
plot(ts_analytical, mass)
%
mass(end) / mass_stationary(end)
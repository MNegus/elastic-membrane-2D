%% validation_plotting.m
% Conducts validation by comparing the FD solution to the exact solution
% for varying grid sizes and timestep sizes

close all;

% Parent directory where all of the data is stored
parent_dir = "forced_validation_data";

%% N_MEMBRANE validation
MAX_TIMESTEP = 2499; % Max value of the timestep in all the dirs

% Values of N_MEMBRANE
N_MEMBRANES = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096];

N_errors = zeros(size(N_MEMBRANES));

for k = 1 : length(N_MEMBRANES)
   N_MEMBRANE = N_MEMBRANES(k);
   
   data_dir = sprintf("%s/N_MEMBRANE_%d", parent_dir, N_MEMBRANE);
   
   max_error = 0;
   
   for q = 0 : MAX_TIMESTEP
      sprintf("N_MEMBRANE = %d, q = %d", N_MEMBRANE, q)
      implicit_mat = dlmread(sprintf("%s/mitchell_outputs/w_%d.txt", data_dir, q));  
      exact_mat = dlmread(sprintf("%s/exact_outputs/w_%d.txt", data_dir, q));  
      
      max_diff = max(abs(exact_mat(:, 2) - implicit_mat(:, 2)))
      if max_error < max_diff
          max_error = max_diff;
      end
   end
   
   N_errors(k) = max_error;
end

%% Plot N_MEMBRANE data
close(figure(1));
figure(1);
loglog(N_MEMBRANES, N_errors, '-o');
hold on;
plot(N_MEMBRANES, 0.25 ./ N_MEMBRANES.^2);
xlabel("Number of nodes, N");
ylabel("Max norm error");
ylim([0.5 * min(N_errors), 1.5 * max(N_errors)]);
grid on;
legend(["Data", "y ~ 1 / N^2"]);
title("Data for DT = 1e-4");
exportgraphics(gcf, "grid_size_convergence.png");

%% Timestep validation

% Values o N_MEMBRANE
N_MEMBRANE = 4096;

% Timestep powers
DTS = [1e-1 3e-2 1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5];
tmax = 0.25;
timestep_errors = zeros(size(DTS));
max_idxs = zeros(size(timestep_errors));

for k = 1 : length(DTS)
    DT = DTS(k);
    
    MAX_TIMESTEP = floor(tmax / DT - 1);
    
    data_dir = sprintf("%s/DT_%d", parent_dir, k);

    max_error = 0;
    max_idx = 0;

    for q = 0 : MAX_TIMESTEP
%       sprintf("DT = %g, q = %d", DT, q)
%       implicit_mat = dlmread(sprintf("%s/implicit_outputs/w_%d.txt", data_dir, q));
      implicit_mat = dlmread(sprintf("%s/mitchell_outputs/w_%d.txt", data_dir, q)); 
      exact_mat = dlmread(sprintf("%s/exact_outputs/w_%d.txt", data_dir, q));  

      [max_diff, I] = max(abs(exact_mat(:, 2) - implicit_mat(:, 2)));
      if max_error < max_diff
          max_error = max_diff;
          max_idx  = I;
      end
      q
      timestep_errors
    end

%     % Compare just w_1
%     q = 1;
%     implicit_mat = dlmread(sprintf("%s/mitchell_outputs/w_%d.txt", data_dir, q)); 
%     exact_mat = dlmread(sprintf("%s/exact_outputs/w_%d.txt", data_dir, q));  
%     max_error = max(abs(exact_mat(:, 2) - implicit_mat(:, 2)))


    
    timestep_errors(k) = max_error
%     max_idxs(k) = max_idx
end

%% Plot timestep errors
timestep_errors
max_idxs
close(figure(2));
figure(2);
loglog(DTS, timestep_errors, '-o');
hold on;
plot(DTS, 40 * DTS.^2);
set(gca, 'XDir','reverse');
xlabel("Timestep, dt");
ylabel("Max norm error");
title("Data for N = 2048");
grid on;
legend("Data", "y ~ dt^2");
exportgraphics(gcf, "timestep_convergence.png");



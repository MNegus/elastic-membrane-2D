%% condition_number.m
% Code to plot the condition number of the matrix A for varying parameters


%% Parameters
ALPHA = 0.1; BETA = 1; GAMMA = 1; 
L = 4;
N_MEMBRANES = [512, 1024, 2048, 4096, 8192];
legend_entries = strings(size(N_MEMBRANES));
for n = 1 : length(N_MEMBRANES)
    legend_entries(n) = sprintf("N = %d\n", N_MEMBRANES(n));
end
T_MAX = 0.1;
DELTA_POWERS = linspace(-1, -5, 5);
DELTA_TS = 10.^DELTA_POWERS;

%% Loop over dts and dxs 
close(figure(1));
figure(1);
hold on;

for N_MEMBRANE = N_MEMBRANES
    N_MEMBRANE
    M = N_MEMBRANE - 1; % We ignore the end point
    xs = linspace(0, L, N_MEMBRANE);
    Deltax = L / (N_MEMBRANE - 1); 
    
    conds = zeros(size(DELTA_TS));

    for k = 1 : length(DELTA_TS)
        
        % Setting DELTA_T
        DELTA_T = DELTA_TS(k)
        Cbeta = BETA * Deltax^2 / GAMMA;
        Calpha = ALPHA * Deltax^4 / (GAMMA * DELTA_T^2);

        % A definition
        A_upper_upper = ones(M, 1);
        A_upper_upper(3) = 2;

        A_upper = (-Cbeta - 4) * ones(M, 1);
        A_upper(2) = -2 * Cbeta - 8;

        A_main = (4 * Calpha + 2 * Cbeta + 6) * ones(M, 1);
        A_main(2) = A_main(2) + 1;
        A_main(M) = A_main(M) - 1;

        A_lower = (-Cbeta - 4) * ones(M, 1);

        A_lower_lower = ones(M, 1);

        A = spdiags([A_lower_lower A_lower A_main A_upper A_upper_upper], -2:2, M, M);
        
        % Determine cond(A) and save
        conds(k) = condest(A);
    end
    
    plot(DELTA_TS, conds, '-o');
    
end

%%
set(gca, 'XDir','reverse');
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
xlabel("dt");
ylabel("cond(A)");
legend(legend_entries, 'location', 'northeast');
function etas = jet_root_etas(x, no_points)
%% jet_root_etas
% Finds the values of etas in a straight line from x_tilde = x / J(t) up to
% the free surface
% ONLY WORKS FOR VALUES UP TO x = 1

    % Solver options
    options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);

    %% Function definition to find eta_0
    eta0_fun = @(eta) x - (1 / pi) * real(eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
    
    %% Find eta0, the value of eta along the substrate
    eta0_guess = -0.1; % eta0 will be between -1 and 0
    eta0 = fsolve(eta0_fun, eta0_guess, options)
    
    %% Find the value of z, using the equation for the free surface
    eta1_fun = @(eta) x - (1 / pi) * (eta - log(eta) - 1);
    eta1_guess = 0.5; % eta1 will be between 0 and 1
    eta1 = fsolve(eta1_fun, eta1_guess, options)
    height = 1 + 4 * sqrt(eta1) / pi;
    
    
    %% Defines range of etas to search through
    %%
    num_reals = 1e3;
    num_imags = 1e3;
    eta_reals = linspace(eta0 - 0.1, eta1 + 0.1, num_reals);
    eta_imags = linspace(0, 0.5, num_imags);

    %% Finds x and z values
    xs = [];
    zs = [];
    eta_points = [];
    
%     close all;
%     figure(1);
%     eta_heights = linspace(0, 1, 1e3);
%     x_heights = (1 / pi) * (eta_heights - log(eta_heights) - 1);
%     z_heights = 1 + 4 * sqrt(eta_heights) / pi;
    
    for eta_imag = eta_imags
        for eta_real = eta_reals
        
            eta = eta_real + 1i * eta_imag;
            zeta = (1 / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);

            xs(end + 1) = real(zeta);
            zs(end + 1) = imag(zeta);
            eta_points(end + 1) = eta;
        end 
 
    end
%     plot(x_heights, z_heights);
%     hold on
%     scatter(xs, zs);
%     hold off;
%     drawnow;

    %% Find query points in z
    zq = linspace(0, height, no_points);
    xq = x * ones(size(zq));

    etas = griddata(xs, zs, eta_points, xq, zq);

    %% Find z and x values from etas and plot them
%     zetas = (1 / pi) * (etas + 4 * 1i * sqrt(etas) - log(etas) + 1i * pi - 1);
%     hold on;
%     scatter(real(zetas), imag(zetas));
%     drawnow;
    
end
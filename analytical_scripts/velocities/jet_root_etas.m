function etas = jet_root_etas(x, no_points)
%% jet_root_etas
% Finds the values of etas in a straight line from x_tilde = x / J(t) up to
% the free surface

    % Solver options
    options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Steptolerance', 1e-8, 'OptimalityTolerance', 1e-10);

    %% Function definition to find eta_0
    eta_0_fun = @(eta) x - (1 / pi) * real(eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
    
    %% Find eta0, the value of eta along the substrate
    eta0_guess = -0.1; % eta0 will be between -1 and 0
    eta0 = fsolve(eta_0_fun, eta0_guess, options)
    
    %% Find the value of z, using the equation for the free surface
    
    
    
%     %%
%     num_reals = 1e2;
%     num_imags = 1e2;
%     eta_reals = linspace(eta0, 1.1, num_reals);
%     eta_imags = linspace(0, 0.5, num_imags);
% 
%     %% Finds x and z values
%     xs = [];
%     zs = [];
%     eta_points = [];
%     for eta_real = eta_reals
%         for eta_imag = eta_imags
%             eta = eta_real + 1i * eta_imag;
%             zeta = (1 / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);
% 
%             xs(end + 1) = real(zeta);
%             zs(end + 1) = imag(zeta);
%             eta_points(end + 1) = eta;
%         end 
%     end
% 
%     %% Find query points in z
%     zq = linspace(0, 1 + 4 / pi, no_points);
%     xq = zeros(size(zq));
% 
%     etas = griddata(xs, zs, eta_points, xq, zq);
% 
%     
    etas = 0
    
end
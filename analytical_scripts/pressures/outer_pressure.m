 function ps = outer_pressure(xs, m_tt_fun, w_tt_fun, d, A, epsilon)
    ps = zeros(size(xs));
    
    if (d == 0)
        return;
    end
    
    xhats = xs / epsilon;

    %% Finds idx such that xhat = d
    idx = sum(xhats < d);
    
    %% Determines the singular integral for 0 <= x < epsilon * d
%     xhat_dependents = zeros(idx, 1);
    s_vals = xhats(1 : idx)';
    
%     for m = 1 : idx
%         xhat = xhats(m);
%         
%         %% trapz method (USE CUMTRAPZ)
%         integrand =  trapz_integrand(s_vals, xhat, m, d, m_tt_fun, w_tt_fun, epsilon);
%         
%         integral_value = trapz(s_vals, integrand, 2);
%         xhat_dependents(m) = xhat * m_tt_fun(xhat) - integral_value;
% 
%     end
%     
    %% Alt, using only one call to trapz but still using a loop
%     xhat_dependents_alt = zeros(idx, 1); 
%     integrands = zeros(idx, idx);
%     for m = 1 : idx
%         xhat = xhats(m);
%         integrands(m, :) = trapz_integrand(s_vals, xhat, m, d, m_tt_fun, w_tt_fun, epsilon);
%         xhat_dependents_alt(m) = xhat * m_tt_fun(xhat);
%     end
%     integrals = trapz(s_vals, integrands, 2);
%     xhat_dependents_alt = xhat_dependents_alt - integrals;
%     abs(xhat_dependents_alt - xhat_dependents)
%     
    %% Full matrix method
    xhat_dependents_alt = zeros(idx, 1);
    for m = 1 : idx
        xhat = xhats(m);
        xhat_dependents_alt(m) = xhat * m_tt_fun(xhat);
    end
    
    full_integrands = trapz_matrix(s_vals, xhats(1 : idx), d, m_tt_fun, w_tt_fun, epsilon);
    integrals = trapz(s_vals, full_integrands, 2);
    xhat_dependents_alt = xhat_dependents_alt - integrals;
    
%     abs(xhat_dependents_alt - xhat_dependents)
    
    %% Returns p values
%     ps(1 : idx) = (A + xhat_dependents) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    ps(1 : idx) = (A + xhat_dependents_alt) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    
    %% Function definitions
    function vals = trapz_integrand(s, xhat, xhat_idx, d, m_tt_fun, w_tt_fun, epsilon)

        % Sets the bulk nodes
        vals = (2 / pi) * (sqrt(d^2 - s.^2) .* (s .* m_tt_fun(s) - xhat * m_tt_fun(xhat))) ./ (s.^2 - xhat^2);

        % Adjusts the node where s = xhat to avoid singularity
        if (xhat == 0)
            vals(xhat_idx) = (2 / pi) * d * w_tt_fun(0);
        else 
            vals(xhat_idx) = (2 / pi) * sqrt(d^2 - xhat^2) * (m_tt_fun(xhat) + xhat * w_tt_fun(epsilon * xhat)) / (2 * xhat);
        end 
    end

    function mat = trapz_matrix(s_vals, xhats, d, m_tt_fun, w_tt_fun, epsilon)
        
        % Sets the bulk nodes
        mat = (2 / pi) * (sqrt(d^2 - s_vals.^2) .* (s_vals .* m_tt_fun(s_vals) ...
            - xhats .* m_tt_fun(xhats))) ./ (s_vals.^2 - (xhats).^2);
        
        % Set the diagonal to zero (as it would be NaN)
        mat(1 : 1 + size(mat,1) : end) = 0;
        
        % Sets the diagonals
        diagonal = (2 / pi) * sqrt(d^2 - xhats.^2) .* (m_tt_fun(xhats) + xhats .* w_tt_fun(epsilon * xhats)) ./ (2 * xhats);
        mat = mat + diag(diagonal);
        mat(1, 1) = (2 / pi) * d * w_tt_fun(0);
    end

end
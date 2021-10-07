 function ps = outer_pressure_flat(xs, w_t_fun, w_tt_fun, d, d_t, epsilon)
    ps = zeros(size(xs));
    
    idx = sum(epsilon * xs < d);
    xhats = xs(1 : idx) / epsilon;
    
    ps(1 : idx) = (1 / epsilon) * (-w_tt_fun(epsilon * xhats) .* sqrt(d^2 - xhats.^2) ...
        + ((1 - w_t_fun(epsilon * xhats)) * d * d_t) ./ sqrt(d^2 - xhats.^2));
    
 
 end
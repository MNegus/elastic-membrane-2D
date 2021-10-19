

%% Params
EPSILON = 1;
t = 1e-4;
d = 2 * sqrt(t);
d_t = 1 / sqrt(t);
A = d * d_t;
tol = 1e-4;
xs_adapt = adaptive_x_values(d, A, tol, EPSILON)
d_idx = sum(xs_adapt <= d);
ps = zeros(size(xs_adapt));
ps(1 : d_idx) = (1 / EPSILON) * A ./ sqrt(d^2 - xs_adapt(1 : d_idx).^2);



plot(xs_adapt, ps, '-o');
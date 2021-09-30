alpha = 2;
beta = 1;
gamma = 2;
epsilon = 1;
L = 16;
N = 128;

lambda = @(n) pi * (2 * n - 1) / (2 * L);
as = zeros(N, 1);

xs = linspace(0, L, 1024);
dt = 1e-2;
tmax = 0.25;
tvals = 0 : dt : tmax


close(figure(1));
figure(1);
for q = 1 : length(tvals)
   t = tvals(q)
   
   for n = 1 : N
        k = beta * lambda(n)^2 + gamma * lambda(n)^4;
        as(n) = 10 * (-2 * alpha + k * t^2 + 2 * alpha * cos(sqrt(k / alpha) * t)) / (k^2 * sqrt(n));
   end
   
   plot(xs, w_solution(xs, as, L, N));
   drawnow;
   pause(0.1);
    
end


function ws = w_solution(xs, as, L, N)
    ws = zeros(size(xs));
    
    lambda = @(n) pi * (2 * n - 1) / (2 * L);
    
    
    %% FIND AN ALTERNATIVE USING MATRIX MULTIPLICATION
    for n = 1 : N
        ws = ws + as(n) * cos(lambda(n) * xs) / sqrt(L);
    end

end


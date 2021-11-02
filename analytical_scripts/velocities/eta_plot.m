close all;

figure(1);
for x = [0, 0.01, 0.1, 1]
    etas = jet_root_etas(x, 256);
    scatter(real(etas), imag(etas));
    hold on;
end

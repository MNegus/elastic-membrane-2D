clear;
close all;

fontsize = 22;
%% Red-blue colour map
cmap = [   59 76 192;...
        68 90 204;...
        77 104 215;...
        87 117 225;...
        98 130 234;...
        108 142 241;...
        119 154 247;...
        130 165 251;...
        141 176 254;...
        152 185 255;...
        163 194 255;...
        174 201 253;...
        184 208 249;...
        194 213 244;...
        204 217 238;...
        213 219 230;...
        221 221 221;...
        229 216 209;...
        236 211 197;...
        241 204 185;...
        245 196 173;...
        247 187 160;...
        247 177 148;...
        247 166 135;...
        244 154 123;...
        241 141 111;...
        236 127 99;...
        229 112 88;...
        222 96 77;...
        213 80 66;...
        203 62 56;...
        192 40 47;...
        180 4 38] / 255;


%%
DELTA_T = 1e-4;
q = 10;
EPSILON = 1;
T_MAX = 0.4 - 0.125;
L = 16;
k = @(n) pi * (2 * n - 1) / (2 * L);

%% Physical parameters

Ns = [16, 32, 64, 128, 256, 512, 1024, 2048];

color_idxs = floor(linspace(1, length(cmap), length(Ns)));

tiledlayout(1, 2, 'TileSpacing','Compact', 'Padding', 'Compact');

for varying = ["alpha", "betas"]
    nexttile;
    hold on;
    
    if varying == "alpha"
        ALPHAS = 2.^linspace(-7, 4, 1e3);
        BETAS = 0;
        INDEP_PARAM = ALPHAS;
    else
        ALPHAS = 1;
        BETAS = 10 * 2.^linspace(-2, 8, 1e3);
        INDEP_PARAM = BETAS;
    end
    
    %% Areas
    min_gamma = 1e-4;
    for N_idx = 1 : length(Ns)
        N = Ns(N_idx);

        %% Stability solution
        gammas_stable = ALPHAS * (2 * pi / (EPSILON^2 * q * DELTA_T))^2 / k(N)^4 - BETAS / k(N)^2;

        %% Reflection solution
        gammas_pos = (ALPHAS * L^2 - T_MAX^2 * BETAS * EPSILON^4 ...
            + sqrt(L^2 * ALPHAS .* (L^2 * ALPHAS + 2 * T_MAX^2 * BETAS * EPSILON^4)))...
            / (2 * k(N)^2 * T_MAX^2 * EPSILON^4);
        gammas_negs = (ALPHAS * L^2 - T_MAX^2 * BETAS * EPSILON^4 ...
            - sqrt(L^2 * ALPHAS .* (L^2 * ALPHAS + 2 * T_MAX^2 * BETAS * EPSILON^4)))...
            / (2 * k(N)^2 * T_MAX^2 * EPSILON^4); 

        %% Set upper and lower bounds
        min_gamma = 1e-6;
        upper_gammas = min(gammas_pos, gammas_stable);
        lower_gammas = max(gammas_negs, min_gamma * ones(size(gammas_negs)));

        %% Find logical values
        Li = (lower_gammas <= upper_gammas);

        patch([INDEP_PARAM(Li) fliplr(INDEP_PARAM(Li))], [lower_gammas(Li) fliplr(upper_gammas(Li))], cmap(color_idxs(N_idx), :), ...
            'edgecolor', 'none')

        plot(INDEP_PARAM, gammas_pos, 'color', cmap(color_idxs(N_idx), :), ...
            'linewidth', 2, 'linestyle', ':');
        plot(INDEP_PARAM, gammas_stable, 'color', cmap(color_idxs(N_idx), :), ...
            'linewidth', 2, 'linestyle', '--');
        plot(INDEP_PARAM, upper_gammas, 'color', cmap(color_idxs(N_idx), :), ...
            'linewidth', 2);
    end

    %%
    ylim([10^-5, 10^4]);
    yticks(10.^(-4 : 2 :  5));
    xlim([min(INDEP_PARAM), max(INDEP_PARAM)])

    set(gca, 'Yscale', 'log');
    set(gca, 'Xscale', 'log');
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    
    if varying == "alpha"
        xticks(2.^(-6 : 3 : 3));
        xlabel("$\alpha$", "interpreter", "latex", "Fontsize", fontsize);
        ylabel("$\gamma$", "interpreter", "latex", "Fontsize", fontsize);
        
        alpha_scatter = 2.^[0, 0.5, 1, 1.5, 2, 2.5, 3];
        scatter(alpha_scatter, 2 * ones(size(alpha_scatter)), [], 'black');
        
        gamma_scatter = [0.5, 1, 2, 4, 8, 16, 32];
        scatter(2 * ones(size(gamma_scatter)), gamma_scatter, [], '+', 'black');
        
        title("$\beta$ = 0", "interpreter", "latex", "fontsize", fontsize);
    else
        xticks(10 * 2.^(-2 : 3 : 8));
        xlabel("$\beta$", "interpreter", "latex", "Fontsize", fontsize);
        ylabel("$\gamma$", "interpreter", "latex", "Fontsize", fontsize);
%         set(gca,'YTickLabel',[]);
        
        beta_scatter = [5, 10, 20, 40, 80, 160, 320, 640, 1280];
        scatter(beta_scatter, 16 * ones(size(beta_scatter)), [], 'd', 'black');
        
        title("$\alpha$ = 1", "interpreter", "latex", "fontsize", fontsize);
    end
    Ax = gca;
    Ax.Layer = 'top';
    
%     legend();
    
end

x0=400;
y0=400;
width=1000;
height=500;

set(gcf,'position',[x0,y0,width,height])

exportgraphics(gcf, "figures/parameter_contours.png", "Resolution", 300);
savefig(gcf, "figures/parameter_contours.fig");


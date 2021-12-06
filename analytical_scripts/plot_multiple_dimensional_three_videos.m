%% plot_multiple_fd.m

close all;
clear;

%% Dimensional parameters
R = 1e-3;
V = 5;


%% Parameters
[EPSILON, ~, ~, ~, L, T_MAX, DELTA_T, N_MEMBRANE] ...
    = parameters();

no_params = 7;

% Spatial parameters
DELTA_X = L / (N_MEMBRANE - 1); 
xs = (0 : DELTA_X : L - DELTA_X)';

% Basilisk parameters
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = 0;
T_VALS = 0 : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;

% Pressure type (composite or outer)
pressure_type = "composite";


%% Data dirs
parent_dir = "/media/michael/newarre/elastic_membrane/parameter_sweeping";

%% Sets up figure
figure(1);
% layout.Padding = 'compact';
fontsize = 20;

x0=400;
y0=400;
width=1200;
height=249;

set(gcf,'position',[x0,y0,width,height])
set(gcf,'color','w');

%% Set up animated lines
figure(1);
membrane_lines = arrayfun(@(x) animatedline(), 1:no_params);
xlim([0, 4]);
xticks([0 : 0.5 : 6]);
ylim([-0.08, 0.02]);
ylabel("$-w(x, t)$ (mm)", "interpreter", "latex", "Fontsize", fontsize);
xlabel("$x$ (mm)", "interpreter", "latex", "Fontsize", fontsize);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
grid on;


% colormap("jet")
% cmap = colormap;
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
color_idxs = floor(linspace(1, length(cmap), no_params));
for idx = 1 : no_params
    for line = membrane_lines(no_params -idx + 1)
        color_idxs(idx);
        line.Color = cmap(color_idxs(idx), :);
        line.LineWidth = 2;
    end
end

%% Writer object setup


%% Loops over membrane parameters
varying_idx = 0;
for varying = ["alpha", "beta", "gamma"]
    varying_idx = varying_idx + 1;
    
    membrane_obj = VideoWriter(sprintf('%s_varying_membrane.avi', varying));
    membrane_obj.FrameRate = 30;
    open(membrane_obj);
    
    %% Sets the parameters
    if varying == "alpha"
        ALPHAS = [1, 1.5, 2, 3, 4, 6, 8] / EPSILON^2;
        BETAS = zeros(size(ALPHAS)) * EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "beta"
        BETAS = [0, 10, 40, 160, 640, 2560, 10240] * EPSILON^2;
        ALPHAS = ones(size(BETAS)) / EPSILON^2;
        GAMMAS = 2 * (EPSILON^2 * ALPHAS).^3 * EPSILON^2;
    elseif varying == "gamma"
        GAMMAS = [2, 8, 32, 128, 512, 2048, 8192] * EPSILON^2;
        ALPHAS = 2 * ones(size(GAMMAS)) / EPSILON^2;
        BETAS = zeros(size(GAMMAS)) * EPSILON^2;
    end
    
    for k = 1 : 20 : length(T_VALS)
        %% Updates time
        t = T_VALS(k);
        t

        if (t <= 0) 
            continue 
        end
        
        %% Loops over parameters and plots
        
        for idx = 1 : no_params
            ALPHA = ALPHAS(no_params - idx + 1);
            BETA = BETAS(no_params - idx + 1);
            GAMMA = GAMMAS(no_params - idx + 1);

            % Loads in parameters
            parameter_dir = sprintf("%s/%s_varying/alpha_%g-beta_%g-gamma_%g/finite_differences/%s", ...
              parent_dir, varying, ALPHA, BETA, GAMMA, pressure_type);
            displayname = ['$\alpha =$ ', num2str(ALPHA),', $\beta =$ ', num2str(BETA), ', $\gamma =$ ', num2str(GAMMA)];


            ws_mat = matfile(sprintf("%s/w_%d.mat", parameter_dir, k - IMPACT_TIMESTEP));
            ws = EPSILON^2 * ws_mat.w_next;

            % w plot  
            clearpoints(membrane_lines(idx));
            addpoints(membrane_lines(idx), xs, -ws);
            
            title(sprintf("$t$ = %.4f (ms)", t / V), "Interpreter", "latex"); 
   
        end
        
        f1 = getframe(gcf);
        writeVideo(membrane_obj, f1);

        drawnow;
    end
    
    
end
close(membrane_obj);
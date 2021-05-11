elastic-membrane-2D

Overall this repo aims to model 2D droplet impact onto an elastic membrane. The position of the elastic membrane at z = -w(x, t) is given by the PDE

ALPHA w_tt - BETA w_xx + GAMMA w_xxxx = p(x, t),

where p(x, t) is the pressure along the membrane. The repo is split in these directories:

finite_differences: Code to solve the membrane equation using finite difference methods. This contains code (in C or MATLAB) for the algorithms, but for a given pressure function p(x, t) in order to develop and verify the methods
normal_modes: Code to solve the membrane equation using the normal modes method, usually purely in MATLAB. This has the algorithms for given pressure functions, but not fully coupled
droplet_impact: Code to solve the coupled problem where the pressure p(x, t) is taken from the pressure in the droplet responding to the motion of the membrane. This is either an analytical model taken from Wagner theory, or a numerical model using VOF simulations.

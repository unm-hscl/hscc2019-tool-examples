%
% Name        : Figure2.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
% Date        : 2018-10-11
%
% Description : Generate Figure 2 from submitted work; plot 2-d set comparison
%               of stochastic viability analysis for a double integrator to
%               stay within a box [-1, 1]^2 using SReachTools, specifically 
%               SReachSet and SReachDynProg
% 

close all;
clearvars;

% example parameters
T = 0.25;
time_horizon = 5;

NO_OF_DIR_VECS_CC = 24;
NO_OF_DIR_VECS_GP = 24;

% warnings off to supress output
warning('off','all');

% safe set definition
safe_set = Polyhedron('lb', -1 * ones(1, 2), 'ub', ones(1, 2));
% target tube definition
target_tube = Tube('viability', safe_set, time_horizon);

sys = getChainOfIntegLtiSystem(2, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
    RandomVector('Gaussian', zeros(2, 1), diag([1e-6, 1e-3])));

fprintf('    Computation times\n');
fprintf('    -----------------\n');

% Lagrangian
% -------------
fprintf('    Lagrangian Underapproximation: ');
opts = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'box', 'err_thresh', 1e-3);

tic;
luSet = SReachSet('term', 'lag-under', sys, 0.8, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);

fprintf('    Lagrangian Overapproximation: ');
opts = SReachSetOptions('term', 'lag-over', ...
    'bound_set_method', 'box', 'err_thresh', 1e-3);

tic;
loSet = SReachSet('term', 'lag-over', sys, 0.8, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);

% Convex chance-constrained set methods
% -----------------------------------------
theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS_CC);
set_of_direction_vectors = [cos(theta_vector); 
                            sin(theta_vector)];

fprintf('    Convex Chance-Constrained: ');
opts = SReachSetOptions('term', 'chance-open', 'pwa_accuracy', 1e-3, ...
    'set_of_dir_vecs', set_of_direction_vectors,...
    'init_safe_set_affine',Polyhedron());

tic;
cccSet = SReachSet('term', 'chance-open', sys, 0.8, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);

% FT with Genz and PatternSearch
% ----------------------------------------
theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS_GP);
set_of_direction_vectors = [cos(theta_vector); 
                            sin(theta_vector)];

fprintf('    Fourier Transform Genz PatternSearch: ');
opts = SReachSetOptions('term', 'genzps-open', 'desired_accuracy', 1e-3, ...
    'set_of_dir_vecs', set_of_direction_vectors,...
    'init_safe_set_affine',Polyhedron(),'verbose',0);

tic;
genzSet = SReachSet('term', 'genzps-open', sys, 0.8, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);
fprintf('\n');

% Dynamic Programing
% ---------------------
fprintf('    Dynamic Programming: ');
tic;
[prob_x, cell_of_xvec] = SReachDynProg('term', sys, 0.05, 0.01, target_tube);
ct = toc;
fprintf('%.5f\n', ct)
fprintf('        --> x_inc = 0.05\n');
fprintf('        --> u_inc = 0.01\n');
fprintf('\n');

dyn_soln_lvl_set = getDynProgLevelSets2D(cell_of_xvec, prob_x, 0.8, ...
    target_tube);

figure()
plot(safe_set, 'color', [0.95, 0.95, 0]);
hold on;
plot(loSet, 'color', 'r');
plot(dyn_soln_lvl_set, 'color', 'b');
plot(cccSet, 'color', [1, 0.6, 0]);
plot(genzSet, 'color', [0, 0.6, 1]);
plot(luSet, 'color', 'g','alpha',0.5);
leg = legend('Safe set', 'lag-over', 'SReachDynProg', 'chance-open', ...
    'genzps-open', 'lag-under');
set(leg,'Location','EastOutside');
hold off;
axis square;
box on;
grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
title('Stochastic reach sets computed at $\alpha=0.8$','interpreter','latex');

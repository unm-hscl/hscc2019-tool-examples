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
prob_thresh = 0.8;

NO_OF_DIR_VECS_CC = 24;
NO_OF_DIR_VECS_GP = 24;

% warnings off to supress output
% warning('off','all');

% safe set definition
safe_set = Polyhedron('lb', -1 * ones(1, 2), 'ub', ones(1, 2));
% target tube definition
target_tube = Tube('viability', safe_set, time_horizon);

sys = getChainOfIntegLtiSystem(2, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
    RandomVector('Gaussian', zeros(2, 1), diag([1e-6, 1e-3])));

fprintf('    Computation times\n');
fprintf('    -----------------\n');

%% Lagrangian
% -------------
fprintf('    Lagrangian Underapproximation: ');
opts = SReachSetOptions('term', 'lag-under', 'bound_set_method','ellipsoid',...
   'compute_style','vfmethod','vf_enum_method','lrs');

tic;
luSet = SReachSet('term', 'lag-under', sys, prob_thresh, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);

fprintf('    Lagrangian Overapproximation: ');
opts = SReachSetOptions('term', 'lag-over', 'bound_set_method','ellipsoid',...
   'compute_style','vfmethod','vf_enum_method','lrs');

tic;
loSet = SReachSet('term', 'lag-over', sys, prob_thresh, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);

%% Convex chance-constrained set methods
% -----------------------------------------
theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS_CC);
set_of_direction_vectors = [cos(theta_vector); 
                           sin(theta_vector)];

fprintf('    Convex Chance-Constrained: ');
opts = SReachSetOptions('term', 'chance-open', 'pwa_accuracy', 1e-3, ...
   'set_of_dir_vecs', set_of_direction_vectors,...
   'init_safe_set_affine',Polyhedron());

tic;
cccSet = SReachSet('term', 'chance-open', sys, prob_thresh, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);

%% FT with Genz and PatternSearch
% ----------------------------------------
theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS_GP);
set_of_direction_vectors = [cos(theta_vector); 
                           sin(theta_vector)];

fprintf('    Fourier Transform Genz PatternSearch: ');
opts = SReachSetOptions('term', 'genzps-open', 'desired_accuracy', 5e-2, ...
   'set_of_dir_vecs', set_of_direction_vectors, ... 'PSoptions', optimset('Display','iter'), 
   'init_safe_set_affine',Polyhedron(),'verbose', 1, 'tol_bisect', 1e-3);


tic;
genzSet = SReachSet('term', 'genzps-open', sys, prob_thresh, target_tube, opts);
ct = toc;
fprintf('%.5f\n', ct);
fprintf('\n');

%% Dynamic Programing
% ---------------------
x_step = 0.05;
u_step = 0.01;
fprintf('    Dynamic Programming: ');
tic;
[prob_x, cell_of_xvec] =SReachDynProg('term', sys, x_step, u_step, target_tube);
ct = toc;
fprintf('%.5f\n', ct)
fprintf('        --> x_inc = %1.2f\n', x_step);
fprintf('        --> u_inc = %1.2f\n', u_step);
fprintf('\n');

dyn_soln_lvl_set = getDynProgLevelSets2D(cell_of_xvec, prob_x, prob_thresh, ...
    target_tube);

%% Plot all the sets
figure(1)
clf;
plot(safe_set, 'color', [0.95, 0.95, 0]);
hold on;
plot(loSet, 'color', 'r');
% plot(dyn_soln_lvl_set, 'color', 'b');
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
title(sprintf('Stochastic reach sets computed at $\\alpha=%1.2f$', ...
    prob_thresh),'interpreter','latex');

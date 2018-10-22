%
% Name        : chainOfIntsExample.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
% Date        : 2018-10-13
%
% Description : Generate Figure 3 from submitted work; scalability plot of
%               different SReachSet methods on a chain of integrators system.
% 
% Notes
%   - For Lagrangian methods we only show scalability comoparison of the 
%     underapproximation methods
% 

clearvars;
close all;

SCALABILITY_MAT_NAME = 'scalability_comptimes.mat';
LAG_SCALE_LIMIT = 7;
CCC_SCALE_LIMIT = 20;
GENZPS_SCALE_LIMIT = 20;
NO_OF_DIR_VECS = 24;

% -----------------------------------------
% 
% Global parameters
% 
% ------------------------------------------

% example parameters
T = 0.25;
time_horizon = 5;

% turn off warning to supress output
warning('off', 'all');

% -----------------------------------------
% 
% Lagrangian Methods
% 
% -----------------------------------------

fprintf('Chain of Integrators: Lagrangian approximations\n');
fprintf('-----------------------------------------------\n\n');


lag_comp_times = zeros(1, LAG_SCALE_LIMIT - 1);
for lv = 2:LAG_SCALE_LIMIT
    fprintf('    Dimension: %d\n', lv);

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(lv, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(lv,1), 0.001*eye(lv)));

    opts = SReachSetOptions('term', 'lag-under', ...
        'bound_set_method', 'box', 'err_thresh', 1e-3);
    % S = SReachSetLagBset(sys.dist, 0.95, opts);

    tic;
    luSet = SReachSet('term', 'lag-under', sys, 0.8, target_tube, opts);
    lag_comp_times(lv - 1) = toc;

    fprintf('    Computation Time: %.5f\n', lag_comp_times(lv-1));
    fprintf('\n');
end

lag = struct('comptimes', lag_comp_times, 'run_time', datestr(now));
save(SCALABILITY_MAT_NAME, 'lag', '-append');

% -------------------------------------------------
% 
% Convex chance-constrained methods
% 
% -------------------------------------------------

fprintf('Chain of Integrators: Convex chance-constrained approximations\n');
fprintf('-----------------------------------------------\n\n');

ccc_comp_times = zeros(1, CCC_SCALE_LIMIT - 1);
for lv = 2:CCC_SCALE_LIMIT
    fprintf('    Dimension: %d\n', lv);

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    mu = zeros(lv, 1);
    sig = diag([1e-6*ones(1, lv-1), 1e-3]);

    sys = getChainOfIntegLtiSystem(lv, T, ...
        Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', mu, sig));

    theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS);
    set_of_direction_vectors = [cos(theta_vector); 
                                sin(theta_vector);
                                zeros(lv-2, length(theta_vector))];

    init_safe_set_affine = Polyhedron('He',...
        [zeros(lv-2,2) eye(lv-2) zeros(lv-2,1)]);
    opts = SReachSetOptions('term', 'chance-open', 'pwa_accuracy', 1e-3, ...
        'set_of_dir_vecs', set_of_direction_vectors,...
        'init_safe_set_affine', init_safe_set_affine);

    tic;
    cccSet = SReachSet('term', 'chance-open', sys, 0.8, target_tube, opts);
    ccc_comp_times(lv-1) = toc;

    fprintf('    Computation Time: %.5f\n', ccc_comp_times(lv-1));
    fprintf('\n');
end

ccc = struct('comptimes', ccc_comp_times, 'run_time', datestr(now));
save(SCALABILITY_MAT_NAME, 'ccc', '-append');


% -----------------------------------------------------
% 
% Fourier Transform with Genz Algorithm
% 
% -----------------------------------------------------

fprintf('Chain of Integrators: Genz-ps approximations\n');
fprintf('-----------------------------------------------\n\n');

genzps_comp_times = zeros(1, GENZPS_SCALE_LIMIT - 1);
for lv = 2:GENZPS_SCALE_LIMIT
    fprintf('    Dimension: %d\n', lv);

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    mu = zeros(lv, 1);
    sig = diag([1e-6*ones(1, lv-1), 1e-3]);

    sys = getChainOfIntegLtiSystem(lv, T, ...
        Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', mu, sig));

    theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS);
    set_of_direction_vectors = [cos(theta_vector); 
                                sin(theta_vector);
                                zeros(lv-2, length(theta_vector))];

    init_safe_set_affine = Polyhedron('He',...
        [zeros(lv-2,2) eye(lv-2) zeros(lv-2,1)]);
    opts = SReachSetOptions('term', 'genzps-open',...
        'desired_accuracy', 1e-3, ...
        'set_of_dir_vecs', set_of_direction_vectors(:,1:1:end),...
        'init_safe_set_affine', init_safe_set_affine,'verbose',1);

    tic;
    genzpsSet = SReachSet('term','genzps-open',sys, 0.8, target_tube, opts);
    genzps_comp_times(lv-1) = toc;

    fprintf('    Computation Time: %.5f\n', genzps_comp_times(lv-1));
    fprintf('\n');
end

genzps = struct('comptimes', genzps_comp_times, 'run_time', datestr(now));
save(SCALABILITY_MAT_NAME, 'genzps', '-append');


hf = figure(3);
plot(2:length(lag_comp_times)+1, lag_comp_times, ...
    'Color', 'b', ...
    'Marker', '^', ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'b');
hold on;
plot(2:length(ccc_comp_times)+1, ccc_comp_times, ...
    'Color', 'k', ...
    'Marker', 'x', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k');
plot(2:length(genzps_comp_times)+1, genzps_comp_times, ...
    'Color', 'r', ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r');
hold off;

ha = gca;
ha.YScale = 'log';
ha.YLim = [10e-3, 10e4];

lh = legend('lag-under', 'chance-open', 'genzps-open');
lh.Location = 'southeast';

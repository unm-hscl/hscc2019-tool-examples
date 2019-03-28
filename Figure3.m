%
% Name        : Figure3.m
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
LAG_O_SCALE_LIMIT = 12;
DP_SCALE_LIMIT = 4;
LAG_U_SCALE_LIMIT = 5;
CCC_SCALE_LIMIT = 20;
GENZPS_SCALE_LIMIT = 20;
CC_NO_OF_DIR_VECS = 24;
GP_NO_OF_DIR_VECS = 24;
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
%% Lag-under
fprintf('Chain of Integrators: Lagrangian (under)approximations\n');
fprintf('------------------------------------------------------\n\n');

lag_u_comp_times = zeros(1, LAG_U_SCALE_LIMIT - 1);
for lv = 2:LAG_U_SCALE_LIMIT
    fprintf('    Dimension: %d\n', lv);

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(lv, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(lv,1), 0.001*eye(lv)));
    n_dim = sys.state_dim + sys.input_dim;
    opts = SReachSetOptions('term', 'lag-under', 'bound_set_method',...
        'ellipsoid','compute_style', 'support','system', sys, ...
        'n_vertices', 2^n_dim * 6 + 2*n_dim, 'vf_enum_method', 'lrs');

    tic;
    luSet = SReachSet('term', 'lag-under', sys, 0.8, target_tube, opts);
    lag_u_comp_times(lv - 1) = toc;

    fprintf('    Computation Time: %.5f\n', lag_u_comp_times(lv-1));
    fprintf('\n');
end

lagunder = struct('comptimes', lag_u_comp_times, 'run_time', datestr(now));
% save(SCALABILITY_MAT_NAME, 'lagunder', '-append');
save(SCALABILITY_MAT_NAME, 'lagunder');

%% Lag-over
fprintf('Chain of Integrators: Lagrangian (over)approximations\n');
fprintf('-----------------------------------------------------\n\n');

lag_o_comp_times = zeros(1, LAG_O_SCALE_LIMIT - 1);
for lv = 2:LAG_O_SCALE_LIMIT
    fprintf('    Dimension: %d\n', lv);

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(lv, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(lv,1), 0.001*eye(lv)));

    n_dim = sys.state_dim;
    opts = SReachSetOptions('term', 'lag-over', 'bound_set_method',...
        'ellipsoid','compute_style','support','system', sys, ...
        'n_vertices', 2^n_dim * 6 + 2*n_dim, 'vf_enum_method', 'lrs');

    tic;
    SReachSet('term', 'lag-over', sys, 0.8, target_tube, opts);
    lag_o_comp_times(lv - 1) = toc;

    fprintf('    Computation Time: %.5f\n', lag_o_comp_times(lv-1));
    fprintf('\n');
end

lagover = struct('comptimes', lag_o_comp_times, 'run_time', datestr(now));
save(SCALABILITY_MAT_NAME, 'lagover', '-append');

% -------------------------------------------------
% 
% Convex chance-constrained methods
% 
% -------------------------------------------------

fprintf('Chain of Integrators: Convex chance-constrained approximations\n');
fprintf('--------------------------------------------------------------\n\n');

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

    theta_vector = linspace(0, 2*pi, CC_NO_OF_DIR_VECS);
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

%% Fourier Transform
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

    theta_vector = linspace(0, 2*pi, GP_NO_OF_DIR_VECS);
    set_of_direction_vectors = [cos(theta_vector); 
                                sin(theta_vector);
                                zeros(lv-2, length(theta_vector))];

    init_safe_set_affine = Polyhedron('He',...
        [zeros(lv-2,2) eye(lv-2) zeros(lv-2,1)]);
    opts = SReachSetOptions('term', 'genzps-open',...
        'desired_accuracy', 5e-2, ...
        'set_of_dir_vecs', set_of_direction_vectors(:,1:1:end),...
        'init_safe_set_affine', init_safe_set_affine,'verbose', 0);

    tic;
    genzpsSet = SReachSet('term','genzps-open',sys, 0.8, target_tube, opts);
    genzps_comp_times(lv-1) = toc;

    fprintf('    Computation Time: %.5f\n', genzps_comp_times(lv-1));
    fprintf('\n');
end

genzps = struct('comptimes', genzps_comp_times, 'run_time', datestr(now));
save(SCALABILITY_MAT_NAME, 'genzps', '-append');

% -----------------------------------------
% 
% Dynamic programming
% 
% -----------------------------------------
%% Dynamic programming
fprintf('Dynamic programming\n');
fprintf('-----------------\n\n');

dp_comp_times = zeros(1, DP_SCALE_LIMIT - 1);
% Step sizes for gridding
dyn_prog_xinc = 0.1; % For testing, use 0.5
dyn_prog_uinc = 0.05; % For testing, use 0.1

for lv = 2:DP_SCALE_LIMIT
    fprintf('    Dimension: %d\n', lv);
    
    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(lv, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(lv,1), 0.001*eye(lv)));

    tic;
    [prob_x, cell_of_xvec] =  SReachDynProg('term', sys, dyn_prog_xinc, ...
        dyn_prog_uinc, target_tube);
    dp_comp_times(lv - 1) = toc;

    fprintf('    Computation Time: %.5f\n', dp_comp_times(lv-1));
    fprintf('\n');
end

dp = struct('comptimes', dp_comp_times, 'run_time', datestr(now));
save('./MatFiles/scalability_comptimes_02132019.mat', 'dp', '-append')

%% Plotting
hf = figure(3);
plot(2:length(lag_u_comp_times)+1, lag_u_comp_times, ...
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

sendToMatlabTrigger('scalability-hscc.txt', ...
    sprintf('Simulation is done at %s'), datestr(now));

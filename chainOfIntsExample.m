%
% Name        : chainOfIntsExample.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
%
% Description : Script to execute relavant code to the chain of integrators
%     example in the HSCC 2019 Tool paper about SReachTools
% 

clearvars;
close all;

SCALABILITY_MAT_NAME = 'scalability_comptimes.mat';
LAG_SCALE_LIMIT = 7;
CCC_SCALE_LIMIT = 20;
GENZPS_SCALE_LIMIT = 5;
NO_OF_DIR_VECS = 24;

run_methods = {'set'};%'ccc','genzps'};

run_lag = false;
run_dp = false;
run_ccc_set = false;
run_genzps_open_set = false;
run_2d_set = false;

if ischar(run_methods) && strcmpi(run_methods, 'all')
    run_lag = true;
    run_dp = true;
    run_ccc_set = true;
    run_genzps_open_set = true;
    run_2d_set = true;
else
    if iscell(run_methods)
        for lv = 1:length(run_methods)
            switch(lower(run_methods{lv}))
                case 'dynamic-programming'
                    run_dp = true;
                case 'dp'
                    run_dp = true;
                case 'lagrangian'
                    run_lag = true;
                case 'ccc'
                    run_ccc_set = true;
                case 'genzps'
                    run_genzps_open_set = true;
                case 'set'
                    run_2d_set = true;
                otherwise
                    if strcmpi(run_methods{lv}, 'all')
                        error(['When specifying to run all methods, ', ...
                            'run_methods should be a character array ''all''']);
                    else
                        error('Unhandled run_method option');
                    end
            end
        end
    end
end

% -----------------------------------------
% 
% Global parameters
% 
% ------------------------------------------

% example parameters
T = 0.25;
time_horizon = 5;

% -------------------------------------------
% 
% 2d Set Plot (All Methods)
% 
% -------------------------------------------

if run_2d_set
    % Plot 2d double integrator set comparison methods

    warning('off','all');

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, 2), 'ub', ones(1, 2));
    % target tube definition
    target_tube = TargetTube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(2, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(2, 1), diag([1e-6, 1e-3])));

    % Lagrangian
    % -------------
    disp('lag-under');
    opts = SReachSetOptions('term', 'lag-under', ...
        'bound_set_method', 'box', 'err_thresh', 1e-3);

    luSet = SReachSet('term', 'lag-under', sys, 0.8, target_tube, opts);

    disp('lag-over');
    opts = SReachSetOptions('term', 'lag-over', ...
        'bound_set_method', 'box', 'err_thresh', 1e-3);

    loSet = SReachSet('term', 'lag-over', sys, 0.8, target_tube, opts);

    % Dynamic Programing
    % ---------------------
    disp('dynprog');
    [prob_x, cell_of_xvec] = SReachDynProg('term', sys, 0.05, 0.1, target_tube);

    dyn_soln_lvl_set = getDynProgLevelSets2D(cell_of_xvec, prob_x, 0.8, ...
        target_tube);

    % Convex chance-constrained set methods
    % -----------------------------------------
    theta_vector = linspace(0, 2*pi, NO_OF_DIR_VECS);
    set_of_direction_vectors = [cos(theta_vector); 
                                sin(theta_vector)];

    disp('chance-open');
    opts = SReachSetOptions('term', 'chance-open', 'pwa_accuracy', 1e-3, ...
        'set_of_dir_vecs', set_of_direction_vectors,...
        'init_safe_set_affine',Polyhedron());

    cccSet = SReachSet('term', 'chance-open', sys, 0.8, target_tube, opts);

    % FT with Genz and PatternSearch
    % ----------------------------------------
    disp('genzps-open');
    opts = SReachSetOptions('term', 'genzps-open', 'desired_accuracy', 1e-3, ...
        'set_of_dir_vecs', set_of_direction_vectors,...
        'init_safe_set_affine',Polyhedron(),'verbose',0);
    genzSet = SReachSet('term', 'genzps-open', sys, 0.8, target_tube, opts);

    figure()
    plot(safe_set, 'color', [0.95, 0.95, 0]);
    hold on;
    plot(loSet, 'color', 'r');
    plot(dyn_soln_lvl_set, 'color', 'b');
    plot(cccSet, 'color', [1, 0.6, 0]);
    plot(genzSet, 'color', [0, 0.6, 1]);
    plot(luSet, 'color', 'g','alpha',0.5);
    hold off;
    axis square;
    box on;
    grid on;
    xlabel('$x_1$','interpreter','latex');
    ylabel('$x_2$','interpreter','latex');
    warning('on','all');

end


% -----------------------------------------
% 
% Lagrangian Methods
% 
% -----------------------------------------

if run_lag

    % determine underapproximation of sets using lagrangian methods for chain
    % of integrators up to 6-d

    fprintf('Chain of Integrators: Lagrangian approximations\n');
    fprintf('-----------------------------------------------\n\n');

    warning('off','all');

    lag_comp_times = zeros(1, LAG_SCALE_LIMIT - 1);
    for lv = 2:LAG_SCALE_LIMIT
        fprintf('    Dimension: %d\n', lv);

        % safe set definition
        safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
        % target tube definition
        target_tube = TargetTube('viability', safe_set, time_horizon);

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

    warning('on', 'all');

end

% ------------------------------------------------
% 
% Dynamic Programming
% 
% ------------------------------------------------

if run_dp
    % dynamic programming can only really handle chain of integrators of up to
    % 3 dimensions but the simulation time will be long. Going to compute
    % in 2 dimensions for several different grid sizes

    fprintf('Chain of Integrators: Dynamic Programming\n');
    fprintf('-----------------------------------------------\n\n');

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, 2), 'ub', ones(1, 2));
    % target tube definition
    target_tube = TargetTube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(2, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(2,1), 0.001*eye(2)));

    fprintf('    x_inc: 0.05\n');
    fprintf('    u_inc: 0.1\n');
    tic;
    [prob_x, cell_of_xvec] = SReachDynProg('term', sys, 0.02, 0.1, target_tube);
    ct = toc;
    fprintf('    Computation time: %.5f\n\n', ct);

    dyn_soln_lvl_set = getDynProgLevelSets2D(cell_of_xvec, prob_x, 0.8, ...
        target_tube);

    figure()
    plot(safe_set, 'color', [0.8, 0.8, 0]);
    hold on;
    plot(dyn_soln_lvl_set, 'color', 'b');
    hold off;
end

% -------------------------------------------------
% 
% Convex chance-constrained methods
% 
% -------------------------------------------------

if run_ccc_set
    % Convex chance-constrained optimization

    fprintf('Chain of Integrators: Convex chance-constrained approximations\n');
    fprintf('-----------------------------------------------\n\n');

    warning('off','all');

    ccc_comp_times = zeros(1, CCC_SCALE_LIMIT - 1);
    for lv = 2:CCC_SCALE_LIMIT
        fprintf('    Dimension: %d\n', lv);

        % safe set definition
        safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
        % target tube definition
        target_tube = TargetTube('viability', safe_set, time_horizon);

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

    warning('on', 'all');
end

% -----------------------------------------------------
% 
% Fourier Transform with Genz Algorithm
% 
% -----------------------------------------------------

if run_genzps_open_set
    fprintf('Chain of Integrators: Genz-ps approximations\n');
    fprintf('-----------------------------------------------\n\n');

    warning('off','all');

    genzps_comp_times = zeros(1, GENZPS_SCALE_LIMIT - 1);
    for lv = 2:GENZPS_SCALE_LIMIT
        fprintf('    Dimension: %d\n', lv);

        % safe set definition
        safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
        % target tube definition
        target_tube = TargetTube('viability', safe_set, time_horizon);

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
            'init_safe_set_affine', init_safe_set_affine,'verbose',0);

        tic;
        genzpsSet = SReachSet('term','genzps-open',sys, 0.8, target_tube, opts);
        genzps_comp_times(lv-1) = toc;

        fprintf('    Computation Time: %.5f\n', genzps_comp_times(lv-1));
        fprintf('\n');
    end

    genzps = struct('comptimes', genzps_comp_times, 'run_time', datestr(now));
    save(SCALABILITY_MAT_NAME, 'genzps', '-append');

    warning('on', 'all');
end

% figure();
% plot(safe_set, 'color', 'k', 'alpha',1);
% hold on;
% % plot(loSet, 'color', 'y');
% plot(luSet, 'color', 'g');
% hold off;
% xlabel('$x_1$', 'Interpreter', 'latex')
% ylabel('$x_2$', 'Interpreter', 'latex')
% box on;
% leg=legend('Safe set','Overapproximation','Underapproximation');
% set(leg,'Location','EastOutside');
%
% Name        : cwhExample.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
%
% Description : Script to execute relavant code to the chain of integrators
%     example in the HSCC 2019 Tool paper about SReachTools
% 

clearvars;
close all;

run_methods = {'lagrangian'};

run_lag = false;
run_dp = false;
run_ccc_set = false;
run_genzps_open_set = false;

if ischar(run_methods) && strcmpi(run_methods, 'all')
    run_lag = true;
    run_dp = true;
    run_ccc_set = true;
    run_genzps_open_set = true;
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
                case 'genz'
                    run_genzps_open_set = true;
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
time_horizon = 3;                                             % Stay within a line of sight cone for 4 time steps and 
                                                            % reach the target at t=5% Safe Set --- LoS cone
%% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and |vy|<=vymax
zmax = 5;
vmax = 0.5 * ones(1, 3);
A_safe_set = [ 1, 1,  0,  0,  0,  0;           
              -1, 1,  0,  0,  0,  0; 
               0, 0,  1,  0,  0,  0;
               0, 0, -1,  0,  0,  0;
               0, 0,  0,  1,  0,  0;
               0, 0,  0, -1,  0,  0;
               0, 0,  0,  0,  1,  0;
               0, 0,  0,  0, -1,  0;
               0, 0,  0,  0,  0,  1;
               0, 0,  0,  0,  0, -1];
b_safe_set = [0;
              0;
              zmax;
              zmax;
              vmax(1);
              vmax(1);
              vmax(2);
              vmax(2);
              vmax(3);
              vmax(3);];

safe_set = Polyhedron(A_safe_set, b_safe_set);
%% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
target_set = Polyhedron('lb', [-0.1; -0.1; -0.1; -0.01; -0.01; -0.01],...
                        'ub', [0.1; 0; 0.1; 0.01; 0.01; 0.01]);
%% Target tube
target_tube = TargetTube('reach-avoid', safe_set, target_set, time_horizon);

umax = 0.1;
mean_disturbance = zeros(6,1);
covariance_disturbance = diag([1e-4, 1e-4, 1e-4, 5e-8, 5e-8, 5e-8]);


sys = getCwhLtiSystem(6,...
    Polyhedron('lb', -umax*ones(3,1), 'ub',  umax*ones(3,1)), ...
    RandomVector('Gaussian', mean_disturbance, covariance_disturbance));


% -----------------------------------------
% 
% Lagrangian Methods
% 
% -----------------------------------------

if run_lag
    % determine underapproximation of sets using lagrangian methods
    % Because the CWH system is separable we can split its dynamics into a
    % 4-d (x, y, \dot{x}, \dot{y}) and 2-d (z, \dot{z}) system to improve 
    % computation. Otherwise the 6-d system is at the edge of Lagrangian 
    % methods ability to compute because of the vertex-fact problem.

    sys4d = LtiSystem('StateMatrix', [], ...
        'InputMatrix', [], ...
        'InputSpace', Polyhedron('lb', -umax*ones(2,1), 'ub',  umax*ones(2,1))
        'DisturbanceMatrix', eye(4), ...
        'Disturbance', RandomVector('Gaussian', zeros(4, 1), ...
            diag([1e-4, 1e-4, 5e-8, 5e-8])));

    fprintf('Spacecraft Rendezvous Docking: Lagrangian approximations\n');
    fprintf('----------------------------------------------------------\n\n');

    warning('off','all');

    opts = SReachSetOptions('term', 'lag-under', ...
            'bound_set_method', 'box', 'err_thresh', 1e-3);

    tic;
    luSet = SReachSet('term', 'lag-under', sys, 0.8, target_tube, opts);
    ct = toc;

    fprintf('    Computation Time: %.5f\n', ct);
    fprintf('\n');

    warning('on', 'all');

end

% ------------------------------------------------
% 
% Dynamic Programming
% 
% ------------------------------------------------

if run_dp
    % The CWH system is 6-dimensional and is outside of the realm of dynamic
    % programming
end

% -------------------------------------------------
% 
% Convex chance-constrained methods
% 
% -------------------------------------------------

if run_ccc_set
        % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, 2), 'ub', ones(1, 2));
    % target tube definition
    target_tube = TargetTube('viability', safe_set, time_horizon);

    sys = getChainOfIntegLtiSystem(2, T, Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(2,1), 0.001*eye(2)));

    % Convex chance-constrained set methods
    opts = SReachSetOptions('term', 'chance-open', 'pwa_accuracy', 1e-3);

    cccSet = SReachSet('term', 'chance-open', sys, 0.8, target_tube, opts);

end

% -----------------------------------------------------
% 
% Fourier Transform with Genz Algorithm
% 
% -----------------------------------------------------

if run_genzps_open_set

end

figure()
plot(safe_set.slice([3, 4, 5, 6], [0, 0, 0, 0]), 'color', [0.8, 0.8, 0]);
hold on;
plot(target_set.slice([3, 4, 5, 6], [0, 0, 0, 0]), 'color', [0.8, 0.8, 0]);
plot(luSet.slice([3, 4, 5, 6], [0, 0, 0, 0]), 'color', 'b');
hold off;
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
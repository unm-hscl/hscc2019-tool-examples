%
% Name        : Figure5.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
% Date        : 2018-10-12
%
% Description : Generate Figure 5 from submitted work; verification of satellite
%               rendezvous-docking problem using Clohessy-Wiltshire-Hill
%               dynamics
% 

close all;
clearvars;

%% System definition
umax = 0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the
% chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4, Polyhedron('lb', -umax*ones(2,1),...
                                    'ub',  umax*ones(2,1)),...
       RandomVector('Gaussian', mean_disturbance,covariance_disturbance));

%% Target tube construction --- reach-avoid specification
time_horizon=5;          % Stay within a line of sight cone for 4 time steps and 
                         % reach the target at t=5% Safe Set --- LoS cone
% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and 
% |vy|<=vymax
ymax = 2;
vxmax = 0.5;
vymax = 0.5;
A_safe_set = [1, 1, 0, 0;           
             -1, 1, 0, 0; 
              0, -1, 0, 0;
              0, 0, 1,0;
              0, 0,-1,0;
              0, 0, 0,1;
              0, 0, 0,-1];
b_safe_set = [0;
              0;
              ymax;
              vxmax;
              vxmax;
              vymax;
              vymax];
safe_set = Polyhedron(A_safe_set, b_safe_set);

% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01],...
                        'ub', [0.1; 0; 0.01; 0.01]);
target_tube = Tube('reach-avoid',safe_set, target_set, time_horizon);                    
slice_at_vx_vy = zeros(2,1);

%% Preparation for set computation
prob_thresh = 0.8;

n_dir_vecs = 16;
theta_vec = linspace(0, 2*pi, n_dir_vecs);
set_of_dir_vecs_ft = [cos(theta_vec);
                      sin(theta_vec);
                      zeros(2,n_dir_vecs)];
n_dir_vecs = 40;
theta_vec = linspace(0, 2*pi, n_dir_vecs);
set_of_dir_vecs_cc_open = [cos(theta_vec);
                           sin(theta_vec);
                           zeros(2,n_dir_vecs)];
init_safe_set_affine = Polyhedron('He',[zeros(2,2) eye(2,2) slice_at_vx_vy]);

%% CC (Linear program approach)
options = SReachSetOptions('term', 'chance-open',...
    'set_of_dir_vecs', set_of_dir_vecs_cc_open,...
    'init_safe_set_affine', init_safe_set_affine);
timer_cc_open = tic;
[polytope_cc_open, extra_info] = SReachSet('term','chance-open', sys,...
    prob_thresh, target_tube, options);  
elapsed_time_cc_open = toc(timer_cc_open);

%% Fourier transform (Genz's algorithm and MATLAB's patternsearch)
options = SReachSetOptions('term', 'genzps-open',...
    'desired_accuracy', 5e-2, 'set_of_dir_vecs', set_of_dir_vecs_ft,...
    'init_safe_set_affine', init_safe_set_affine, 'verbose', 0);
timer_ft = tic;
polytope_ft = SReachSet('term','genzps-open', sys, prob_thresh,...
    target_tube, options);  
elapsed_time_ft = toc(timer_ft);

%% Lagrangian underapproximation
n_dim = sys.state_dim + sys.input_dim;
timer_lagunder_options = tic;
lagunder_options = SReachSetOptions('term', 'lag-under',...
    'bound_set_method', 'ellipsoid', 'compute_style','support',...
    'system', sys, 'vf_enum_method', 'lrs', 'verbose', 0,...
    'n_vertices', 2^n_dim * 6 + 2*n_dim);
elapsed_time_lagunder_options = toc(timer_lagunder_options);

timer_lagunder = tic;
[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

%% Preparation for Monte-Carlo simulations of the optimal controllers
% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e5;
n_sims_to_plot = 5;

%% Plotting and Monte-Carlo simulation-based validation
figure(5);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', [0.95, 0.95, 0]);
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', [0, 0, 0]);
legend_cell = {'Safe set','Target set'};
plot(polytope_cc_open.slice([3,4], slice_at_vx_vy), 'color',[1, 0.6, 0],'alpha', 1);
legend_cell{end+1} = 'Underapprox. polytope (chance-open)';
plot(polytope_ft.slice([3,4], slice_at_vx_vy), 'color',[0, 0.6, 1],'alpha',1);
legend_cell{end+1} = 'Underapprox. polytope (genzps-open)';
plot(polytope_lagunder.slice([3,4], slice_at_vx_vy), 'color',[0.3, 1, 0.3],'alpha',1);
legend_cell{end+1} = 'Underapprox. polytope (lag-under)';
direction_index_to_plot = 38;
if ~isEmptySet(polytope_cc_open)
    init_state = extra_info(2).vertices_underapprox_polytope(:,direction_index_to_plot);
    input_vec = extra_info(2).opt_input_vec_at_vertices(:,direction_index_to_plot);
    opt_reach_avoid = extra_info(2).opt_reach_prob_i(direction_index_to_plot);

    concat_state_realization = generateMonteCarloSims(...
            n_mcarlo_sims,...
            sys,...
            init_state,...
            time_horizon,...
            input_vec);        
    
    % Check if the location is within the target_set or not
    mcarlo_result = target_tube.contains(concat_state_realization);
    [legend_cell] = plotMonteCarlo(' (chance-open)', mcarlo_result,...
        concat_state_realization, n_mcarlo_sims, n_sims_to_plot,...
        sys.state_dim,init_state, legend_cell);
end
legend(legend_cell, 'Location','South');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
box on;
grid on;

%% Compute time
fprintf(['Elapsed time: (genzps-open) %1.3f | (chance-open) %1.3f ',...
    '| (lag-under) %1.3f seconds\n'], elapsed_time_ft, elapsed_time_cc_open, ...
    elapsed_time_lagunder);

%% Helper functions
% Plotting function
function [legend_cell] = plotMonteCarlo(method_str, mcarlo_result,...
    concat_state_realization, n_mcarlo_sims, n_sims_to_plot, state_dim,...
    initial_state, legend_cell)
% Plots a selection of Monte-Carlo simulations on top of the plot

    green_legend_updated = 0;
    red_legend_updated = 0;
    traj_indices = floor(n_mcarlo_sims*rand(1,n_sims_to_plot));
    for realization_index = traj_indices
        % Check if the trajectory satisfies the reach-avoid objective
        if mcarlo_result(realization_index)
            % Assign green triangle as the marker
            plotOptions = {'Color', 'g', 'Marker', '^', ...
                'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g' 'MarkerSize', 5};
            markerString = 'g^';
        else
            % Assign red asterisk as the marker
            plotOptions = {'Color', 'r', 'Marker', 'x', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r' 'MarkerSize', 5};
            markerString = 'rx';
        end

        % Create [x(t_1) x(t_2)... x(t_N)]
        reshaped_X_vector = reshape(...
            concat_state_realization(:,realization_index), state_dim,[]);

        % This realization is to be plotted
        h = plot([initial_state(1), reshaped_X_vector(1,:)], ...
                 [initial_state(2), reshaped_X_vector(2,:)], ...
                 plotOptions{:});
    end
end

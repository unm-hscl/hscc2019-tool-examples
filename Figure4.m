%
% Name        : Figure4.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
% Date        : 2018-10-11
%
% Description : Generate Figure 4 from submitted work; Dubin's vehicle with a 
%               fixed turning sequence, LTV system
% 

close all;
clearvars;

time_horizon = 50;
time_const = 1/2*time_horizon;
init_heading = pi/10;
sampling_time = 0.1;
box_halflength = 4;
omega = pi/time_horizon/sampling_time;
turning_rate = omega*ones(time_horizon,1);
dist_cov = 0.0001;
prob_thresh = 0.9;
no_of_direction_vectors_ccc = 16;
v_nominal = 10;
umax = v_nominal/3*2;
n_mcarlo_sims = 1e3;

%% LTV system definition
[sys, heading_vec] = getDubinsCarLtv('add-dist',...
    turning_rate,...
    init_heading,...
    sampling_time,...
    Polyhedron('lb',0,'ub',umax),...
    eye(2),...
    RandomVector('Gaussian',zeros(2,1), dist_cov * eye(2)));

target_tube_cell = cell(time_horizon + 1,1);

%% Target tube definition
figure(4);
clf;
hold on;
angle_at_the_center = (heading_vec) - pi/2;
center_box = zeros(2, time_horizon + 1);        
for itt=0:time_horizon
    center_box(:, itt+1) = v_nominal * [cos(angle_at_the_center(itt+1))-cos(angle_at_the_center(1));
                                        sin(angle_at_the_center(itt+1))-sin(angle_at_the_center(1))];
    target_tube_cell{itt+1} = Polyhedron('lb',center_box(:, itt+1) - box_halflength * exp(- itt/time_const), 'ub', center_box(:, itt+1) + box_halflength*exp(- itt/time_const));
    plot(target_tube_cell{itt+1},'alpha',0.5,'color','y');
end
axis equal
axis([-8    10   -5   21]);
box on;
grid on;

target_tube = Tube(target_tube_cell{:});

init_state_ccc_open = [2;2] + [-1;1];
init_state_genzps_open = [2;2] + [1;-1];
init_state_particle_open = [2;2] + [0;1];
init_state_ccc_affine = [2;2] + [2;1];

%% Quantities needed to compute the optimal mean trajectory 
% Compute the H matrix
[Z,H,G] = sys.getConcatMats(time_horizon);
% Compute the zero input mean trajectory
sys_no_input = LtvSystem('StateMatrix',sys.state_mat,...
    'DisturbanceMatrix', sys.dist_mat,'Disturbance',sys.dist);
X_zizs = SReachFwd('concat-stoch', sys_no_input,...
    zeros(sys.state_dim,1), time_horizon);
mean_X_zizs = X_zizs.mean();
mean_X_zizs = mean_X_zizs(sys.state_dim+1:end);

%% SReachPoint: chance-open (with pwa_accuracy 1e-3)
opts = SReachPointOptions('term', 'chance-open','pwa_accuracy',1e-3);
tic;
[prob_ccc_open, opt_input_vec_ccc_open] = SReachPoint('term',...
    'chance-open', sys, init_state_ccc_open, target_tube, opts);
elapsed_time_ccc_open = toc;
optimal_mean_X_ccc_open = Z * init_state_ccc_open +...
    H * opt_input_vec_ccc_open + mean_X_zizs;
optimal_mean_trajectory_ccc_open = reshape(optimal_mean_X_ccc_open,...
    sys.state_dim,[]);

%% SReachPoint: chance-affine (with \Delta_u = 0.01)
opts = SReachPointOptions('term', 'chance-affine','max_input_viol_prob',1e-2,...
    'verbose',1);
tic
[prob_ccc_affine, opt_input_vec_ccc_affine, opt_input_gain_ccc_affine] =...
    SReachPoint('term', 'chance-affine', sys, init_state_ccc_affine,...
        target_tube, opts);
elapsed_time_ccc_affine = toc;
% X = Z * x_0 + H * (M \mu_W + d) + G * \mu_W
muW = kron(ones(time_horizon,1), sys.dist.parameters.mean);
optimal_mean_X_ccc_affine = Z * init_state_ccc_affine +...
    H * (opt_input_gain_ccc_affine * muW + opt_input_vec_ccc_affine) + G * muW;
optimal_mean_trajectory_ccc_affine = reshape(optimal_mean_X_ccc_affine,...
    sys.state_dim,[]);

%% SReachPoint: genzps-open
opts = SReachPointOptions('term', 'genzps-open',...
    'PSoptions',psoptimset('display','iter'), 'desired_accuracy', 5e-2);
tic
[prob_genzps_open, opt_input_vec_genzps_open] = SReachPoint('term',...
    'genzps-open', sys, init_state_genzps_open, target_tube, opts);
elapsed_time_genzps = toc;
optimal_mean_X_genzps_open =  Z * init_state_genzps_open +...
    H * opt_input_vec_genzps_open + mean_X_zizs;
optimal_mean_trajectory_genzps_open = reshape(optimal_mean_X_genzps_open,...
    sys.state_dim,[]);

%% SReachPoint: particle-open (use verbosity 1)
tic
opts = SReachPointOptions('term','particle-open','verbose',1,'n_particles',50);
[prob_particle_open, opt_input_vec_particle_open] = SReachPoint('term',...
    'particle-open', sys, init_state_particle_open, target_tube, opts);
elapsed_time_particle = toc;
optimal_mean_X_particle_open =  Z * init_state_particle_open +...
    H * opt_input_vec_particle_open + mean_X_zizs;
optimal_mean_trajectory_particle_open = reshape(optimal_mean_X_particle_open,...
    sys.state_dim,[]);

%% Lagrangian under
timer_lagunder = tic;
theta_polytope_vec = linspace(0,2*pi,10)';
lagunder_options = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'ellipsoid', ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 1);

[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

%% Compute a far-away safe initial
cvx_begin quiet
    variable initial_state(sys.state_dim, 1)
    minimize ([1 1]*initial_state)
    subject to
        polytope_lagunder.A*initial_state <= polytope_lagunder.b;
        target_tube(1).A*initial_state <= target_tube(1).b;
cvx_end
switch cvx_status
    case 'Solved'
        fprintf('Testing initial state: ');
        disp(initial_state');
        
        % Create a controller based on the underapproximation
        srlcontrol = SReachLagController(sys, ... 
            extra_info_under.bounded_dist_set, ...
            extra_info_under.stoch_reach_tube);
        % Generate Monte-Carlo simulations using the srlcontrol and
        % generateMonteCarloSims
        timer_mcarlo = tic;
        [X,U,W] = generateMonteCarloSims(n_mcarlo_sims, sys, ...
            initial_state, time_horizon, srlcontrol, [], ...
            lagunder_options.verbose);
        elapsed_time_mcarlo = toc(timer_mcarlo);
        avg_time_mc = elapsed_time_mcarlo / n_mcarlo_sims;

        % % Plot the convex hull of the spread of the points
        polytopesFromMonteCarloSims(X, 4, [1,2], {'color','k','alpha',0});
        a = gca;            
        for tindx = 1:time_horizon-1
            a.Children(tindx).Annotation.LegendInformation.IconDisplayStyle='off';
        end
        a.Children(1).Annotation.LegendInformation.IconDisplayStyle='on';
        a.Children(1).DisplayName = 'Trajectory spread at various time steps';           
        
        % Plot the initial state
        scatter(initial_state(1), initial_state(2), 200, 'ko', 'filled', ...
            'DisplayName','Initial state');
    otherwise        
end

init_state_lag = initial_state;
optimal_mean_trajectory_lag = reshape(sum(X, 2) / size(X, 2), 2, []);

%% Save data
save_mat_file_path = strcat('./MatFiles/','DubinsCar_example_point_',datestr(now,'YYYYmmDD_HHMMSS'),'.mat');
save(save_mat_file_path);

%% Plot the set
figure(4);
clf;
hold on;
for itt=0:time_horizon
    if itt==0
        % Remember the first the tube
        h_target_tube=plot(target_tube_cell{1},'alpha',0.5,'color','y');
    else
        plot(target_tube_cell{itt+1},'alpha',0.08,'LineStyle',':','color','y');
    end            
end
axis equal        
% Plot the optimal mean trajectory from the vertex under study
h_opt_mean_ccc = scatter(...
      [init_state_ccc_open(1), optimal_mean_trajectory_ccc_open(1,:)],...
      [init_state_ccc_open(2), optimal_mean_trajectory_ccc_open(2,:)],...
      30, 'ko', 'filled','DisplayName', 'chance-open');
h_opt_mean_ccc_affine = scatter(...
      [init_state_ccc_affine(1), optimal_mean_trajectory_ccc_affine(1,:)],...
      [init_state_ccc_affine(2), optimal_mean_trajectory_ccc_affine(2,:)],...
      30, 'gs', 'filled','DisplayName', 'chance-affine');
h_opt_mean_genzps = scatter(...
      [init_state_genzps_open(1), optimal_mean_trajectory_genzps_open(1,:)],...
      [init_state_genzps_open(2), optimal_mean_trajectory_genzps_open(2,:)],...
      30, 'bd', 'filled', 'DisplayName', 'genzps-open');
h_opt_mean_particle = scatter(...
      [init_state_particle_open(1), optimal_mean_trajectory_particle_open(1,:)],...
      [init_state_particle_open(2), optimal_mean_trajectory_particle_open(2,:)],...
      30, 'r^', 'filled','DisplayName', 'particle-open');  
h_opt_mean_lag = scatter(...
      [init_state_lag(1), optimal_mean_trajectory_lag(1,:)],...
      [init_state_lag(2), optimal_mean_trajectory_lag(2,:)],...
      30, 'p', 'filled', 'DisplayName', 'particle-open');  
xlabel('x');
ylabel('y');
axis equal
box on;
legend_cell = {'Target tube', ...
    'chance-open','chance-affine',...
    'genzps-open','particle-open'};
h_vec = [h_target_tube, h_opt_mean_ccc,...
    h_opt_mean_ccc_affine, h_opt_mean_genzps, h_opt_mean_particle];
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');

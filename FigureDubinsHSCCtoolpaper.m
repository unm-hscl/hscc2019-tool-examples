close all
clc
clear

date_str_mat = '20181016_152216';
load(strcat('./MatFiles/DubinsCar_example_',date_str_mat,'.mat'));

% Get a random point
% init_state_ccc_open = ccc_polytope.randomPoint();
init_state_ccc_open = [1.5087;1.7064];
% init_state_genzps_open = ccc_polytope.randomPoint();
init_state_ccc_affine = ccc_polytope.randomPoint();

%% Quantities needed to compute the optimal mean trajectory 
% Compute the H matrix
[Z,H,G] = sys.getConcatMats(time_horizon);
% Compute the zero input mean trajectory
sys_no_input = LtvSystem('StateMatrix',sys.state_mat,...
    'DisturbanceMatrix', sys.dist_mat,'Disturbance',sys.dist);
[mean_X_zizs, ~] = SReachFwd('concat-stoch', sys_no_input,...
    zeros(sys.state_dim,1), time_horizon);

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

% %% SReachPoint: genzps-open (use default options)
% opts = SReachPointOptions('term', 'genzps-open',...
%     'PSoptions',psoptimset('display','iter'));
% [prob_genzps_open, opt_input_vec_genzps_open] = SReachPoint('term',...
%     'genzps-open', sys, init_state_genzps_open, target_tube, opts);
% optimal_mean_X_ccc_open =  Z * init_state_genzps_open +...
%     H * opt_input_vec_genzps_open + mean_X_zizs;
% optimal_mean_trajectory_genzps_open = reshape(optimal_mean_X_genzps_open,...
%     sys.state_dim,[]);

%% SReachPoint: chance-open (with pwa_accuracy 1e-3)
opts = SReachPointOptions('term', 'chance-affine','max_input_viol_prob',1e-2,...
    'verbose',2);
tic
[prob_ccc_affine, opt_input_vec_ccc_affine, opt_input_gain_ccc_affine] =...
    SReachPoint('term', 'chance-affine', sys, init_state_ccc_affine,...
        target_tube, opts);
elapsed_time_ccc_affine = toc;
optimal_mean_X_ccc_affine = Z * init_state_ccc_affine +...
    H * opt_input_gain_ccc_affine *...
        repmat(sys.dist.parameters.mean,time_horizon, 1) +...
    H * opt_input_vec_ccc_open + mean_X_zizs;
optimal_mean_trajectory_ccc_affine = reshape(optimal_mean_X_ccc_affine,...
    sys.state_dim,[]);

%% Plot the set
figure(101);
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
h_nominal_traj=scatter(center_box(1,:), center_box(2,:), 50,'ks','filled');        
% Plot the optimal mean trajectory from the vertex under study
h_opt_mean_ccc = scatter(...
      [init_state_ccc_open(1), optimal_mean_trajectory_ccc_open(1,:)],...
      [init_state_ccc_open(2), optimal_mean_trajectory_ccc_open(2,:)],...
      30, 'bo', 'filled','DisplayName', 'Mean trajectory (chance-open)');
% h_opt_mean_genzps = scatter(...
%       [init_state_genzps_open(1), optimal_mean_trajectory_genzps_open(1,:)],...
%       [init_state_genzps_open(2), optimal_mean_trajectory_genzps_open(2,:)],...
%       30, 'bo', 'filled','DisplayName', 'Mean trajectory (genzps-open)');
h_opt_mean_ccc_affine = scatter(...
      [init_state_ccc_open(1), optimal_mean_trajectory_ccc_affine(1,:)],...
      [init_state_ccc_open(2), optimal_mean_trajectory_ccc_affine(2,:)],...
      30, 'ms', 'filled','DisplayName', 'Mean trajectory (chance-affine)');
xlabel('x');
ylabel('y');
axis equal
axis(axis_vec);
box on;
set(gca,'FontSize',fontSize);
legend_cell = {'Target tube', 'Nominal trajectory',...
    'Mean trajectory (chance-open)','Mean trajectory (chance-affine)'};
h_vec = [h_target_tube, h_nominal_traj, h_opt_mean_ccc, h_opt_mean_ccc_affine];
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');



% %% Check if the location is within the target_set or not
% n_mcarlo_sims = 1e5;
% concat_state_realization = generateMonteCarloSims(n_mcarlo_sims,...
%     sys, init_state_ccc_open, time_horizon, opt_input_vec);
% mcarlo_result = target_tube.contains([repmat(init_state_ccc_open,1,n_mcarlo_sims);
%                                       concat_state_realization]);
% fprintf('SReachPoint prob: %1.2f, Simulated prob: %1.2f', prob, sum(mcarlo_result)/n_mcarlo_sims);
                                  

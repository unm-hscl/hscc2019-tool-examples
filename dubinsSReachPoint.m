dubinsSetup

init_state_ccc_open = [2;2] + [-1;1];

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
xlabel('x');
ylabel('y');
axis equal
axis([-8    10   -5   21]);
box on;
legend_cell = {'Target tube', 'Nominal trajectory', 'Mean trajectory (chance-open)'};
h_vec = [h_target_tube, h_nominal_traj, h_opt_mean_ccc];
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');
set(gca,'FontSize',20);
saveas(gcf, '../results/Dubins_opt_traj.png')

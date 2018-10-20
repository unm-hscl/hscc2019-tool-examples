clear
close all
clc

time_horizon = 50;
time_const = 1/2*time_horizon;
init_heading = pi/10;
sampling_time = 0.1;
box_halflength = 4;
omega = pi/time_horizon/sampling_time;
turning_rate = omega*ones(time_horizon,1);
dist_cov = 0.001;
probability_threshold_of_interest = 0.8;
no_of_direction_vectors_ccc = 16;
v_nominal = 10;
umax = v_nominal/3*2;

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
figure(100);clf;hold on
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

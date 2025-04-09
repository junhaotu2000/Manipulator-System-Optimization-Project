clear; close all hidden;
run setup.m
import casadi.*  

opti_fun = make_arm_opti();
NB = 3;
T = 2.0;   % Total duration (2 seconds)
N = 100;   % Number of time steps
DT = T / N;
time_grid = linspace(0, T, N);
% Define each link's rotational inertia (using diagonal matrices)
rotInertia = [[0.02, 0.02, 0.01]; [0.02, 0.02, 0.01]; [0.02, 0.02, 0.01]];

% Define link masses (vector with 3 elements)
mass = [1.0, 1.0, 1.0];
model = create3DoFRobotModel(rotInertia, mass);

tic
[x, u] = opti_fun(rotInertia, mass);
toc

X_opt = full(x);
U_opt = full(u);

Q_opt  = X_opt(1:NB, :);       % Joint angles (3 x N)
dQ_opt = X_opt(NB+1:end, :);   % Joint velocities (3 x N)

% --- (a) Plot Joint States and Inputs ---
figure('Name','3DOF Robot - States & Inputs','NumberTitle','off');

subplot(1,3,1); % Joint Angles
p1 = plot(time_grid, Q_opt(1,:), 'LineWidth',1.5); hold on;
p2 = plot(time_grid, Q_opt(2,:), 'LineWidth',1.5);
p3 = plot(time_grid, Q_opt(3,:), 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
title('Joint Angles'); grid on;

subplot(1,3,2); % Joint Velocities
plot(time_grid, dQ_opt(1,:), 'LineWidth',1.5); hold on;
plot(time_grid, dQ_opt(2,:), 'LineWidth',1.5);
plot(time_grid, dQ_opt(3,:), 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Velocity (rad/s)');
title('Joint Velocities'); grid on;

subplot(1,3,3); % Joint Torques
plot(time_grid, U_opt(1,:), 'LineWidth',1.5); hold on;
plot(time_grid, U_opt(2,:), 'LineWidth',1.5);
plot(time_grid, U_opt(3,:), 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Torque (Nm)');
title('Joint Torques'); grid on;

hL = axes('Position',[0 0 1 1],'Visible','off');
lgd = legend(hL, [p1, p2, p3], {'Joint1', 'Joint2', 'Joint3'}, ...
    'Orientation','horizontal', 'Location','northoutside');
lgd.Position = [0.35 0.95 0.3 0.05];

% --- (b) Animation of Robot Motion ---
figure('Name','3DOF Robot Animation','NumberTitle','off');
showmotion(model, time_grid, Q_opt);
drawnow;

% Compute the metrics for each joint
avg_torque   = mean(abs(U_opt), 2);         % Average torque [Nm]
peak_velocity = max(abs(dQ_opt), [], 2);      % Peak velocity [rad/s]
final_angle  = Q_opt(:, end);                % Final joint angle [rad]

% Create and display a table with the results
results = table(avg_torque, peak_velocity, final_angle, ...
    'VariableNames', {'Average_Torque_Nm', 'Peak_Velocity_rad_s', 'Final_Angle_rad'}, ...
    'RowNames', {'Joint1', 'Joint2', 'Joint3'});

disp(results);
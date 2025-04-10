clear
close all hidden
run setup.m; % Load the required libraries and functions
import casadi.*

%% 1. Run Mechanical Design Optimization
[rotInertia_arr, mass_out, L] = MechOptimization();
        
fprintf('Mechanical Optimization Complete.\n');
disp('Rotational Inertia Array (each column is [I_x; I_y; I_z] in SI units):');
disp(rotInertia_arr);
disp('Mass Vector (link masses in kg):');
disp(mass_out);

% In the robot planning system the function create3DoFRobotModel expects a 
% 3x3 matrix with each row corresponding to a link's inertia parameters.
% Here the GA returns a 3x3 matrix where each column corresponds to a link.
% Therefore, transpose the array to get the desired format.

% inertias = rotInertia_arr'

%% 2. Run Robot Planning Optimization via CasADi
% L = mechOut.L; % Link lengths
opti_fun = make_arm_opti(); % Create the casadi function for the optimizer
[x, u] = opti_fun(rotInertia_arr', mass_out, L/1000);

X_opt = full(x);
U_opt = full(u);
NB = 3;

% Extract trajectories: joint angles and joint velocities.
Q_opt  = X_opt(1:NB, :);       % Joint angles (3 x N)
dQ_opt = X_opt(NB+1:end, :);    % Joint velocities (3 x N)

fprintf('Robot Planning Optimization Complete.\n');

% Plot the Robot Planning Results
time_grid = linspace(0, 2, size(X_opt,2)); % Assuming 2 sec total duration, N time steps.

figure('Name','Robot States & Inputs','NumberTitle','off');
subplot(1,3,1); % Joint Angles
plot(time_grid, Q_opt(1,:), 'LineWidth', 1.5); hold on;
plot(time_grid, Q_opt(2,:), 'LineWidth', 1.5);
plot(time_grid, Q_opt(3,:), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)'); title('Joint Angles'); grid on;

subplot(1,3,2); % Joint Velocities
plot(time_grid, dQ_opt(1,:), 'LineWidth', 1.5); hold on;
plot(time_grid, dQ_opt(2,:), 'LineWidth', 1.5);
plot(time_grid, dQ_opt(3,:), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (rad/s)'); title('Joint Velocities'); grid on;

subplot(1,3,3); % Joint Torques
plot(time_grid, U_opt(1,:), 'LineWidth', 1.5); hold on;
plot(time_grid, U_opt(2,:), 'LineWidth', 1.5);
plot(time_grid, U_opt(3,:), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Torque (Nm)'); title('Joint Torques'); grid on;

% our rotational inertia is a fixed matrix
% example: rotInertia = [[0.02, 0.02, 0.01]; [0.02, 0.02, 0.01]; [0.02, 0.02, 0.01]];

%% 3. Run Controller Optimization
%The controller optimization is for a single joint, so it will need to run three times.

I_const = 0.02.*ones(length(U_opt(1,:)));

%Joint 1
[RMSE_Final_1, q_1] = run_pid_optimization(I_const, U_opt(1,:), Q_opt(1,:));

%Joint 2
[RMSE_Final_2, q_2] = run_pid_optimization(I_const, U_opt(2,:), Q_opt(2,:));

%Joint 3
[RMSE_Final_3, q_3] = run_pid_optimization(I_const, U_opt(3,:), Q_opt(3,:));

RMSE_Final = [RMSE_Final_1; RMSE_Final_2; RMSE_Final_3]; %RMSE Tracking error 3xN array [joint 1, joint 2, joint 3]
q = [q_1; q_2; q_3]; %Angular position 3xN array [joint 1, joint 2, joint 3]
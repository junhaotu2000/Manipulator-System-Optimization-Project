%% 
% Enhanced Manipulator System Optimization using ATC (Analytical Target Cascading)
% With Multiple Runs and Detailed Result Tables
%

clear all;
close all;
clc;
import casadi.*

q_init = [0;  0;  0]; 
q_goal = [pi/2;  pi/4;  0];
% global opti_fun


opti_fun = make_arm_opti(q_init, q_goal);

%% Top-level system targets
% System-level targets for coupling variables
target_RMSE = 0.05;            % Tracking error target
target_effort = 100;           % Effort minimization target
target_smoothness = 10;        % Trajectory smoothness target
target_path_penalty = 1;       % Path deviation penalty target
target_phase_margin = 30;      % Phase margin target (degrees)

% Combined target vector
system_targets = [target_RMSE; target_effort; target_smoothness; target_path_penalty; target_phase_margin];

% Display system targets
disp('=========== System-level Targets ===========');
disp(['Tracking Error (RMSE): ', num2str(target_RMSE)]);
disp(['Control Effort: ', num2str(target_effort)]);
disp(['Trajectory Smoothness: ', num2str(target_smoothness)]);
disp(['Path Penalty: ', num2str(target_path_penalty)]);
disp(['Phase Margin: ', num2str(target_phase_margin), ' degrees']);
disp('===========================================');


% [X, U, J] = opti_fun(<3x3 matrix of limb inertias around coms>, <3x1 vector of masses>, <3x1 vector of lengths>);
% X = [q1, q2, q3, dq1, dq2, dq3]';
% U = [u1, u2, u3];
% J = f_x;

%% Number of runs to perform
num_runs = 5;

% Initialize result arrays
all_x_star = zeros(11, num_runs);
all_f_star = zeros(1, num_runs);
all_phi_star = zeros(1, num_runs);
all_exitflags = zeros(1, num_runs);
all_mech_vars = zeros(9, num_runs);
all_controller_gains = zeros(3, num_runs);
all_subsystem_metrics = zeros(5, num_runs);

% Timer for performance tracking
tic;

%% Run optimization multiple times
for run_idx = 1:num_runs
    fprintf('\n=========== Starting Run %d of %d ===========\n', run_idx, num_runs);
    
    % Run ATC optimization
    [x_star, f_star, phi_star, exitflag] = manipulator_system_atc_opt(system_targets);
    
    % Store system-level results
    all_x_star(:, run_idx) = x_star;
    all_f_star(run_idx) = f_star;
    all_phi_star(run_idx) = phi_star;
    all_exitflags(run_idx) = exitflag;
    
    % Extract targets for each subsystem
    mech_targets = x_star(1:6);   % Inertia and mass targets
    plan_targets = x_star(7:9);   % Trajectory performance targets
    ctrl_targets = x_star(10:11); % Control performance targets
    
    % Run subsystem optimizations with optimal targets
    fprintf('Solving mechanical subsystem for run %d...\n', run_idx);
    mech_x_star = mechanical_subsystem_atc_opt(mech_targets);
    all_mech_vars(:, run_idx) = mech_x_star;
    
    % Extract link properties from mechanical subsystem solution
    D = mech_x_star(1:3);
    t = mech_x_star(4:6);
    L = mech_x_star(7:9);
    [rotInertia_arr, mass_out, I_trans_arr] = convertLinkParams(D, t, L);
    
    % Simulation parameters
    T = 2.0;            % Total simulation time
    N = 100;            % Number of time steps
    dt = T/N;           % Time step
    q_res = 0.01;       % Position resolution
    tau_max = 15;       % Maximum torque
    
    fprintf('Solving path planning subsystem for run %d...\n', run_idx);
    plan_x_star = planner_subsystem_atc_opt(plan_targets, rotInertia_arr, mass_out, L, opti_fun);
    
    % Extract planned trajectory
    x = reshape(plan_x_star(1:6*N), 6, N);
    q_planned = x(1:3,:);           % Joint position trajectories
    dq_planned = x(4:6,:);          % Joint velocity trajectories
    u_planned = reshape(plan_x_star(6*N+1:end), 3, N);  % Joint torques
    
    % Time-varying inertia and external torque functions
    t_vec = linspace(0, T, N);
    I_func = @(k) 1 + 0.2*sin(t_vec(min(k, N)));  % Time-varying inertia
    tau_ext_func = @(k) 0.1*t_vec(min(k, N));     % Time-varying external torque
    qd_desired = q_planned(1,:);                   % First joint trajectory as reference
    
    fprintf('Solving controller subsystem for run %d...\n', run_idx);
    ctrl_x_star = controller_subsystem_atc_opt(ctrl_targets, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    all_controller_gains(:, run_idx) = ctrl_x_star;
    
    % Calculate final performance metrics
    [RMSE, ~, ~] = pid_simulation(ctrl_x_star, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    effort = sum(u_planned(:).^2);
    smoothness = sum(diff(u_planned,1,2).^2,'all');
    
    % Calculate path penalty
    model = create3DoFRobotModel(rotInertia_arr, mass_out, L/1000);
    path_penalty = 0;
    for k = 1:N-1
        p1 = ee_pos(model, q_planned(:,k));
        p2 = ee_pos(model, q_planned(:,k+1));
        path_penalty = path_penalty + sum((p2-p1).^2);
    end
    
    % Calculate phase margin
    phi_m = phase_margin(ctrl_x_star);
    
    % Store subsystem performance metrics
    all_subsystem_metrics(:, run_idx) = [RMSE; effort; smoothness; path_penalty; phi_m];
    
    % Display condensed run results
    fprintf('Run %d complete. Objective: %.4f, Target deviation: %.4f\n', run_idx, f_star, phi_star);
end

% Record total execution time
total_time = toc;
fprintf('\nAll %d runs completed in %.2f seconds (%.2f minutes)\n', num_runs, total_time, total_time/60);

%% Create and display system-level results table
fprintf('\n=========== System-Level Optimization Results ===========\n');
system_table = table();
system_table.Run = (1:num_runs)';
system_table.Objective = all_f_star';
system_table.TargetDeviation = all_phi_star';
system_table.ExitFlag = all_exitflags';

% Calculate mean and standard deviation
system_table.Properties.UserData.Mean = [mean(all_f_star), mean(all_phi_star)];
system_table.Properties.UserData.Std = [std(all_f_star), std(all_phi_star)];

% Display table
disp(system_table);

% Display statistics
fprintf('Mean Objective: %.4f ± %.4f\n', system_table.Properties.UserData.Mean(1), system_table.Properties.UserData.Std(1));
fprintf('Mean Target Deviation: %.4f ± %.4f\n', system_table.Properties.UserData.Mean(2), system_table.Properties.UserData.Std(2));

%% Create and display system-level design variables table
fprintf('\n=========== System-Level Design Variables ===========\n');
vars_table = table();
vars_table.Run = (1:num_runs)';

% Add each design variable
var_names = {'I_trans_1', 'I_trans_2', 'I_trans_3', 'Mass_1', 'Mass_2', 'Mass_3', ...
             'Effort', 'Smoothness', 'PathPenalty', 'RMSE', 'PhaseMargin'};
         
for i = 1:11
    vars_table.(var_names{i}) = all_x_star(i,:)';
end

% Calculate mean and standard deviation for each variable
means = mean(all_x_star, 2);
stdevs = std(all_x_star, 0, 2);

% Display table
disp(vars_table);

% Display mean and standard deviation
fprintf('\nMean and Standard Deviation of Design Variables:\n');
for i = 1:11
    fprintf('%s: %.6f ± %.6f\n', var_names{i}, means(i), stdevs(i));
end

%% Create and display mechanical subsystem results table
fprintf('\n=========== Mechanical Subsystem Results ===========\n');
mech_table = table();
mech_table.Run = (1:num_runs)';

% Add mechanical design variables
mech_var_names = {'D1', 'D2', 'D3', 't1', 't2', 't3', 'L1', 'L2', 'L3'};
         
for i = 1:9
    mech_table.(mech_var_names{i}) = all_mech_vars(i,:)';
end

% Calculate mean and standard deviation for each variable
mech_means = mean(all_mech_vars, 2);
mech_stdevs = std(all_mech_vars, 0, 2);

% Display table
disp(mech_table);

% Display mean and standard deviation
fprintf('\nMean and Standard Deviation of Mechanical Design Variables:\n');
for i = 1:9
    fprintf('%s: %.2f ± %.2f\n', mech_var_names{i}, mech_means(i), mech_stdevs(i));
end

%% Create and display controller results table
fprintf('\n=========== Controller Subsystem Results ===========\n');
ctrl_table = table();
ctrl_table.Run = (1:num_runs)';

% Add controller gains
ctrl_var_names = {'Kp', 'Ki', 'Kd'};
         
for i = 1:3
    ctrl_table.(ctrl_var_names{i}) = all_controller_gains(i,:)';
end

% Calculate mean and standard deviation for each variable
ctrl_means = mean(all_controller_gains, 2);
ctrl_stdevs = std(all_controller_gains, 0, 2);

% Display table
disp(ctrl_table);

% Display mean and standard deviation
fprintf('\nMean and Standard Deviation of Controller Gains:\n');
for i = 1:3
    fprintf('%s: %.4f ± %.4f\n', ctrl_var_names{i}, ctrl_means(i), ctrl_stdevs(i));
end

%% Create and display performance metrics table
fprintf('\n=========== Performance Metrics Results ===========\n');
metrics_table = table();
metrics_table.Run = (1:num_runs)';

% Add performance metrics
metrics_names = {'RMSE', 'Effort', 'Smoothness', 'PathPenalty', 'PhaseMargin'};
         
for i = 1:5
    metrics_table.(metrics_names{i}) = all_subsystem_metrics(i,:)';
end

% Calculate mean and standard deviation for each metric
metrics_means = mean(all_subsystem_metrics, 2);
metrics_stdevs = std(all_subsystem_metrics, 0, 2);

% Calculate deviation from targets
target_metrics = [target_RMSE; target_effort; target_smoothness; target_path_penalty; target_phase_margin];
metrics_deviation = abs(all_subsystem_metrics - target_metrics);
metrics_deviation_mean = mean(metrics_deviation, 2);

% Display table
disp(metrics_table);

% Display mean and standard deviation
fprintf('\nMean Performance Metrics vs Targets:\n');
for i = 1:5
    fprintf('%s: %.4f ± %.4f (Target: %.4f, Deviation: %.4f)\n', ...
        metrics_names{i}, metrics_means(i), metrics_stdevs(i), target_metrics(i), metrics_deviation_mean(i));
end

%% Visualize the best run results
[~, best_run_idx] = min(all_f_star);
fprintf('\n=========== Visualizing Best Run (Run %d) ===========\n', best_run_idx);

% Extract best run data
best_x_star = all_x_star(:, best_run_idx);
best_mech_vars = all_mech_vars(:, best_run_idx);
best_controller_gains = all_controller_gains(:, best_run_idx);

% Extract link properties
D = best_mech_vars(1:3);
t = best_mech_vars(4:6);
L = best_mech_vars(7:9);
[rotInertia_arr, mass_out, I_trans_arr] = convertLinkParams(D, t, L);

% Simulation parameters
T = 2.0;            % Total simulation time
N = 100;            % Number of time steps
dt = T/N;           % Time step
t_vec = linspace(0, T, N);
q_res = 0.01;       % Position resolution
tau_max = 15;       % Maximum torque
model = create3DoFRobotModel(rotInertia_arr, mass_out, L/1000);

% Plan trajectory
plan_targets = best_x_star(7:9);
plan_x_star = planner_subsystem_atc_opt(plan_targets, rotInertia_arr, mass_out, L,opti_fun);

% Extract trajectory
x = reshape(plan_x_star(1:6*N), 6, N);
q_planned = x(1:3,:);
dq_planned = x(4:6,:);
u_planned = reshape(plan_x_star(6*N+1:end), 3, N);

% Control system
I_func = @(k) 1 + 0.2*sin(t_vec(min(k, N)));
tau_ext_func = @(k) 0.1*t_vec(min(k, N));
qd_desired = q_planned(1,:);
[RMSE, q_actual, error] = pid_simulation(best_controller_gains, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);

% End-effector trajectory
ee_positions = zeros(3, N);
for k = 1:N
    ee_positions(:, k) = ee_pos(model, q_planned(:, k));
end

% Create visualization plots
figure('Name', 'Best Run Results', 'Position', [100, 100, 1200, 800]);

% Joint angles
subplot(2, 3, 1);
plot(t_vec, q_planned', 'LineWidth', 2);
title('Joint Angles');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Joint 1', 'Joint 2', 'Joint 3');
grid on;

% Joint torques
subplot(2, 3, 2);
plot(t_vec, u_planned', 'LineWidth', 2);
title('Joint Torques');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Joint 1', 'Joint 2', 'Joint 3');
grid on;

% End effector path
subplot(2, 3, 3);
plot3(ee_positions(1,:), ee_positions(2,:), ee_positions(3,:), 'b-', 'LineWidth', 2);
hold on;
plot3(ee_positions(1,1), ee_positions(2,1), ee_positions(3,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot3(ee_positions(1,end), ee_positions(2,end), ee_positions(3,end), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('End-Effector Trajectory');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;
legend('Path', 'Start', 'End');

% Control system response
subplot(2, 3, 4);
plot(t_vec, qd_desired, 'b--', t_vec, q_actual, 'r-', 'LineWidth', 2);
title('Control Response (Joint 1)');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Desired', 'Actual');
grid on;

% Tracking error
subplot(2, 3, 5);
plot(t_vec, error, 'k-', 'LineWidth', 2);
title('Tracking Error');
xlabel('Time (s)');
ylabel('Error (rad)');
grid on;
text(0.1, max(error)*0.8, ['RMSE = ', num2str(RMSE)], 'FontSize', 12);

% System metrics
subplot(2, 3, 6);
metrics_data = all_subsystem_metrics(:, best_run_idx);
metrics_names = {'RMSE', 'Effort', 'Smooth', 'Path', 'Phase'};
normalized_metrics = metrics_data ./ [0.1; 500; 100; 10; 60]; % Normalize for visualization
bar(normalized_metrics);
title('Performance Metrics');
xlabel('Metric');
ylabel('Normalized Value');
xticks(1:5);
xticklabels(metrics_names);
grid on;

%% Create a comprehensive summary table for the best run
fprintf('\n=========== Best Run Summary (Run %d) ===========\n', best_run_idx);
summary_table = table();

% System-level results
summary_table.Metric = {'System Objective'; 'Target Deviation'; 'Exit Flag'};
summary_table.Value = [all_f_star(best_run_idx); all_phi_star(best_run_idx); all_exitflags(best_run_idx)];

% System-level design variables
for i = 1:length(var_names)
    summary_table.Metric{end+1} = var_names{i};
    summary_table.Value(end+1) = all_x_star(i, best_run_idx);
end

% Mechanical design variables
for i = 1:length(mech_var_names)
    summary_table.Metric{end+1} = ['Mech_', mech_var_names{i}];
    summary_table.Value(end+1) = all_mech_vars(i, best_run_idx);
end

% Controller gains
for i = 1:length(ctrl_var_names)
    summary_table.Metric{end+1} = ['Ctrl_', ctrl_var_names{i}];
    summary_table.Value(end+1) = all_controller_gains(i, best_run_idx);
end

% Performance metrics
for i = 1:length(metrics_names)
    summary_table.Metric{end+1} = ['Perf_', metrics_names{i}];
    summary_table.Value(end+1) = all_subsystem_metrics(i, best_run_idx);
end

% Display the summary table
disp(summary_table);

%% Save results to MAT file for future analysis
fprintf('\nSaving results to manipulator_atc_results.mat\n');
save('manipulator_atc_results.mat', 'all_x_star', 'all_f_star', 'all_phi_star', ...
     'all_exitflags', 'all_mech_vars', 'all_controller_gains', 'all_subsystem_metrics', ...
     'system_table', 'vars_table', 'mech_table', 'ctrl_table', 'metrics_table', 'summary_table');

%% System-level ATC optimization function
function [x_star, f_star, phi_star, exitflag] = manipulator_system_atc_opt(target)
    % ATC optimization for the complete manipulator system
    % Inputs:
    %   target - 5×1 vector of system-level targets
    % Outputs:
    %   x_star - 11×1 vector of optimal coupling variables
    %   f_star - System-level objective value
    %   phi_star - Target deviation at optimum
    %   exitflag - Optimization exit flag
    
    % Parameters
    weight = [1, 1]; % Weights for [phi, phi_sub_star]
    
    % Design variables: inertia, mass, trajectory, control
    % x = [I_trans_1; I_trans_2; I_trans_3; 
    %      mass_1; mass_2; mass_3;
    %      effort; smoothness; path_penalty;
    %      RMSE; phase_margin]
    
    % Variable bounds
    lb = [0.0001; 0.0001; 0.0001;    % Inertia bounds (kg·m²)
          0.1; 0.1; 0.1;             % Mass bounds (kg)
          0; 0; 0;                   % Trajectory performance bounds
          0; 20];                    % Control performance bounds
      
    ub = [0.1; 0.1; 0.1;             % Inertia bounds (kg·m²)
          10; 10; 10;                % Mass bounds (kg)
          1000; 1000; 100;           % Trajectory performance bounds
          0.2; 60];                  % Control performance bounds
    
    % Initial guess (midpoint)
    x0 = (lb + ub) / 2;
    
    % No linear constraints
    A = []; b = [];
    Aeq = []; beq = [];
    
    % No nonlinear constraints
    nonlcon = [];
    
    % Optimization options
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
                          'MaxFunctionEvaluations', 5000, 'MaxIterations', 50);
    
    % Optimize using ATC helper function
    [x_star, f_star, phi_star, exitflag] = atc_opt(@manipulator_system_atc_obj, target, weight, x0, A, b, Aeq, beq, lb, ub, nonlcon, options, true);
end

%% System-level ATC objective function
function [f, phi, phi_sub_star] = manipulator_system_atc_obj(x, target)
    % ATC objective function for manipulator system
    % Inputs:
    %   x - 11×1 vector of coupling variables
    %   target - 5×1 vector of system-level targets
    % Outputs:
    %   f - System-level objective value
    %   phi - Target deviation
    %   phi_sub_star - Sum of optimal subsystem target deviations
    
    % Extract subsystem targets
    mech_targets = x(1:6);     % Inertia and mass targets
    plan_targets = x(7:9);     % Trajectory performance targets
    ctrl_targets = x(10:11);   % Control performance targets
    
    % Simulation parameters
    dt = 0.02;          % Time step
    T = 2.0;            % Total simulation time
    N = round(T/dt);    % Number of time steps
    q_res = 0.01;       % Position resolution
    tau_max = 15;       % Maximum torque
    
    % Calculate subsystem optimal solutions (keep computation light for system-level)
    [~, ~, mech_phi_star] = mechanical_subsystem_atc_opt(mech_targets);
    
    % Generate placeholder optimization results for planning and control
    % This is a simplified model - only computing deviations not full optimization
    % Get inertia and mass from mechanical targets
    I_trans_targets = mech_targets(1:3);
    mass_targets = mech_targets(4:6);
    
    % Use representative values for L based on mass and inertia relation
    L_approx = zeros(1,3);
    for i = 1:3
        % Approximate L from I and m (using hollow cylinder approximation)
        L_approx(i) = sqrt(12 * I_trans_targets(i) / mass_targets(i)) * 1000; % Convert to mm
    end
    
    % Simplified model for planning and control subsystems
    I_func = @(k) 1 + 0.2*sin(k*dt);
    tau_ext_func = @(k) 0.1*k*dt;
    qd_desired = sin(linspace(0, 4*pi, N));
    
    % Calculate subsystem deviations without full optimizations
    [~, ~, plan_phi_star] = planner_subsystem_simple(plan_targets, diag(I_trans_targets), mass_targets, L_approx);
    [~, ~, ctrl_phi_star] = controller_subsystem_simple(ctrl_targets, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    
    % System-level objective (e.g., total mass + control error)
    f = sum(mass_targets) + ctrl_targets(1);
    
    % Target deviation
    system_performance = [ctrl_targets(1); plan_targets(1); plan_targets(2); plan_targets(3); ctrl_targets(2)];
    phi = sum((target - system_performance).^2);
    
    % Sum of optimal subsystem target deviations
    phi_sub_star = mech_phi_star + plan_phi_star + ctrl_phi_star;
end

%% Simplified planning subsystem for system level calculations
function [x_star, f_star, phi_star] = planner_subsystem_simple(target, rotInertia_arr, mass, L)
    % Compute approximate plan deviation without full optimization
    effort_est = 100 + 10*sum(mass);  % Simplified effort estimate
    smoothness_est = 50;              % Default smoothness
    path_penalty_est = 5;             % Default path penalty
    
    metrics = [effort_est; smoothness_est; path_penalty_est];
    phi_star = sum((target - metrics).^2);
    
    f_star = 0;  % No local objective
    x_star = [];  % No need for decision variables in simplified model
end

%% Simplified controller subsystem for system level calculations
function [x_star, f_star, phi_star] = controller_subsystem_simple(target, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max)
    % Compute approximate controller deviation without full optimization
    K_default = [10; 1; 0.1];  % Default PID gains
    [RMSE, ~, ~] = pid_simulation_quick(K_default, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    phi_m = phase_margin(K_default);
    
    metrics = [RMSE; phi_m];
    phi_star = sum((target - metrics).^2);
    
    f_star = 0;  % No local objective
    x_star = K_default;  % Return default gains
end

%% Quick PID simulation for system level calculations
function [RMSE, q, error] = pid_simulation_quick(K, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max)
    % Simplified PID simulation for fast evaluation
    % Only simulate every 10th step to speed up computation
    step_skip = 10;
    N_reduced = ceil(N/step_skip);
    
    Kp = K(1); Ki = K(2); Kd = K(3);
    q = zeros(1, N);
    dq = zeros(1, N);
    integral_error = 0;
    
    for k = 3:step_skip:N
        e = qd_desired(k) - q(max(1,k-step_skip));
        integral_error = integral_error + e * dt * step_skip;
        derivative_error = e / (dt * step_skip);
        
        tau = Kp * e + Ki * integral_error + Kd * derivative_error;
        tau = max(min(tau, tau_max), -tau_max);
        
        ddq = (tau - tau_ext_func(k)) / I_func(k);
        dq(k) = ddq * dt * step_skip;
        q(k) = q(max(1,k-step_skip)) + dq(k);
        
        % Fill in intermediate steps with linear interpolation
        if k > step_skip
            for j = (k-step_skip+1):(k-1)
                q(j) = q(k-step_skip) + (j-(k-step_skip))*(q(k)-q(k-step_skip))/step_skip;
            end
        end
    end
    
    error = qd_desired - q;
    RMSE = sqrt(mean(error.^2));
end

%% Mechanical Subsystem ATC Optimization (GA-based)
function [x_star, f_star, phi_star, exitflag] = mechanical_subsystem_atc_opt(target)
    % ATC-compatible optimizer for the mechanical subsystem
    % Inputs:
    %   target  - 6×1 vector: [I_trans_target(3 entries); mass_target(3 entries)]
    % Outputs:
    %   x_star  - 9×1 optimal design: [D1;D2;D3; t1;t2;t3; L1;L2;L3]
    %   f_star  - subsystem local objective (unused, = 0)
    %   phi_star- ATC target loss at optimum
    %   exitflag- solver exit flag

    % Get ATC objective function with target as a parameter
    atc_obj = @(x) mechanical_subsystem_atc_obj(x, target);

    % Define the number of design variables (D, t, L)
    nLinks = 3;
    nCont = nLinks * 3;
    
    % Variable bounds for GA (mm)
    lb = [ones(1, nLinks)*20, ones(1, nLinks)*5, ones(1, nLinks)*50];
    ub = [ones(1, nLinks)*200, ones(1, nLinks)*10, ones(1, nLinks)*800];
    
    % Check if we have a fast evaluation (for system-level opt) or full GA (for subsystem)
    if all(target < 1e-10) || any(isnan(target))
        % Fast evaluation for system level
        x_star = (lb + ub)' / 2;  % Middle point for default solution
        [f_star, phi_star, ~] = mechanical_subsystem_atc_obj(x_star, target);
        exitflag = 1;
        return;
    end
    
    % For normal optimization, use pattern search for faster convergence than GA
    % This is more efficient for multiple runs
    options = optimoptions('patternsearch', ...
                          'Display', 'off', ...
                          'MaxFunctionEvaluations', 1000, ...
                          'MeshTolerance', 1e-2);
    
    % Initial point (midpoint)
    x0 = (lb + ub) / 2;
    
    % Run pattern search optimization
    [x_opt, fval_opt, exitflag] = patternsearch(atc_obj, x0, [], [], [], [], lb, ub, @nonlcon_mech, options);
    
    % Format output
    x_star = x_opt';  % Convert to column vector
    [f_star, phi_star, ~] = mechanical_subsystem_atc_obj(x_star, target);
end

%% Mechanical subsystem ATC objective function
function [f, phi, phi_sub_star] = mechanical_subsystem_atc_obj(x, target)
    % Extract design variables (mm)
    nLinks = 3;
    D = x(1:nLinks);
    t = x(nLinks+1:2*nLinks);
    L = x(2*nLinks+1:3*nLinks);
    
    % Fixed carbon fiber material (rho = 1600 kg/m³)
    rho = 1600;
    
    % Calculate actual link masses
    mass_out = zeros(1, nLinks);
    for i = 1:nLinks
        % Unit conversion: mm → m
        D_m = D(i) * 0.001;
        t_m = t(i) * 0.001;
        L_m = L(i) * 0.001;
        
        % Calculate cross-sectional area and volume (hollow circular tube)
        r_outer = D_m / 2;
        r_inner = r_outer - t_m;
        A = pi * (r_outer^2 - r_inner^2);
        V = A * L_m;
        mass_out(i) = rho * V;
    end
    
    % Compute inertia values
    [~, ~, I_trans_arr] = convertLinkParams(D, t, L);
    
    % Local objective (minimize total mass - not used in ATC framework)
    f = sum(mass_out);
    
    % ATC target loss: first 3 entries of target for I_trans_arr, next 3 for mass
    I_trans_targets = target(1:3);
    mass_targets = target(4:6);
    
    phi_inertia = sum((I_trans_targets - I_trans_arr').^2);
    phi_mass = sum((mass_targets - mass_out').^2);
    phi = phi_inertia + phi_mass;
    
    % No deeper subsystems
    phi_sub_star = 0;
end

%% Nonlinear constraints for mechanical design
function [c, ceq] = nonlcon_mech(x)
    % Extract design variables
    nLinks = 3;
    D = x(1:nLinks);                       % Outer diameter (mm)
    t = x(nLinks+1:2*nLinks);              % Wall thickness (mm)
    L = x(2*nLinks+1:3*nLinks);            % Link length (mm)
    
    S = 4;           % Safety factor
    M_max = 300;     % Maximum bending moment (Nm)
    
    % Material properties for carbon fiber (fixed material choice)
    sigma_yield = 500e6;  % Yield strength (Pa)
    E = 60e9;             % Young's modulus (Pa)
    
    c = [];
    for i = 1:nLinks
        % Current link design parameters (units: mm and m)
        D_i = D(i);
        t_i = t(i);
        L_i = L(i);
        
        D_m = D_i * 0.001;
        t_m = t_i * 0.001;
        L_m = L_i * 0.001;
        
        % Calculate cross-sectional inertia
        r_outer = D_m / 2;
        r_inner = r_outer - t_m;
        I = (pi/64) * ((D_m)^4 - (D_m - 2*t_m)^4);
        
        % Calculate bending stress (Pa)
        sigma = (M_max * (D_m/2)) / I;
        % Stress constraint
        c(end+1) = sigma - (sigma_yield / S);
        
        % Manufacturability constraint: ensure the inner diameter is positive (i.e., 2*t_i - D_i <= 0)
        c(end+1) = 2*t_i - D_i;
        
        % Deflection constraint: δ = (M_max * L³) / (3 * E * I) must not exceed 0.1% of L
        delta = (M_max * L_m^3) / (3 * E * I);
        c(end+1) = delta - 0.001 * L_m;
    end
    
    % Workspace constraint: total link length must be at least 1000 mm
    c(end+1) = 1000 - sum(L);
    
    ceq = [];
end

%% Link parameters calculation function
function [rotInertia_arr, mass_out, I_trans_arr] = convertLinkParams(D, t, L, mass_in)
    % Compute rotational inertia and mass for 3 links
    % If mass_in is provided, use it; otherwise calculate from D, t, L
    
    if nargin < 4
        % Calculate mass from link geometry parameters and fixed density
        rho = 1600;  % density kg/m^3 (carbon fiber)
        mass_out = zeros(1, 3);
        
        for i = 1:3
            r_outer = D(i)/2 * 1e-3;        % m
            r_inner = r_outer - t(i)*1e-3;  % m
            
            % Link mass
            volume = pi*(r_outer^2 - r_inner^2) * L(i)*1e-3;  % m^3
            mass_out(i) = rho * volume;                      % kg
        end
    else
        mass_out = mass_in;
    end
    
    % Initialize arrays
    rotInertia_arr = zeros(3, 3);
    I_trans_arr = zeros(1, 3);

    for i = 1:3
        r_outer = D(i)/2 * 1e-3;        % m
        r_inner = r_outer - t(i)*1e-3;  % m
        
        % Rotational inertia
        I_x = 0.5 * mass_out(i) * (r_outer^2 + r_inner^2);
        I_t = (1/12) * mass_out(i) * (3*(r_outer^2 + r_inner^2) + (L(i)*1e-3)^2);
        rotInertia_arr(:, i) = [I_x; I_t; I_t];
        I_trans_arr(i) = I_t;
    end
end

%% Path Planning Subsystem ATC Optimization
function [x_star, f_star, phi_star, exitflag] = planner_subsystem_atc_opt(target, rotInertia_arr, mass_out, L, opti_fun)
    % Inputs:
    %   target      - 3×1 vector of coupling targets: [effort; smoothness; path_penalty]
    %   rotInertia_arr - 3×3 rotational inertia matrix
    %   mass_out    - 1×3 link masses
    %   L           - 1×3 link lengths (mm)
    % Outputs:
    %   x_star      - optimal decision vector
    %   f_star      - local objective
    %   phi_star    - ATC target loss
    %   exitflag    - optimization exit flag
    
    % % Create robot model
    % model = create3DoFRobotModel(rotInertia_arr, mass_out, L/1000);
    % 
    % % Parameters
    % NB = model.NB;
    % N = 100;
    % T = 2.0;
    % DT = T/N;
    % q_init = [0; 0; 0];
    % q_goal = [pi/2; pi/4; 0];
    % tau_max = 15;
    % alpha = 5;
    % 
    % % Initial guess - create a simple trajectory
    % z0 = zeros((2*NB + NB)*N, 1);
    % 
    % % Simple trajectory initialization for better convergence
    % for i = 1:NB
    %     % Create simple trajectory between initial and goal positions
    %     q_traj = linspace(q_init(i), q_goal(i), N);
    % 
    %     % Calculate simple velocity and acceleration
    %     dq_traj = zeros(1, N);
    %     dq_traj(2:end) = diff(q_traj) / DT;
    % 
    %     % Basic torque estimation using inverse dynamics
    %     u_traj = zeros(1, N);
    % 
    %     % Insert into initial guess
    %     for k = 1:N
    %         idx_q = (k-1)*(2*NB) + i;
    %         idx_dq = (k-1)*(2*NB) + NB + i;
    %         idx_u = 2*NB*N + (k-1)*NB + i;
    % 
    %         z0(idx_q) = q_traj(k);
    %         z0(idx_dq) = dq_traj(k);
    %         z0(idx_u) = u_traj(k);
    %     end
    % end
    % 
    % % Optimization options - reduced computation for multi-run
    % options = optimoptions('fmincon', ...
    %     'Algorithm', 'sqp', ...
    %     'Display', 'off', ...
    %     'MaxFunctionEvaluations', 5000, ...
    %     'MaxIterations', 20, ...
    %     'OptimalityTolerance', 1e-3, ...
    %     'StepTolerance', 1e-3, ...
    %     'SpecifyObjectiveGradient', false);
    % 
    % % Create objective function with embedded target
    % cost_func = @(z) planner_subsystem_atc_obj(z, target, NB, N, DT, model, alpha);
    % 
    % % Nonlinear constraints
    % nonlcon_func = @(z) nonlcon_planner(z, NB, N, DT, model, q_init, q_goal, tau_max);
    % 
    % % Run optimization with reduced computation for multi-run
    % [z_opt, fval, exitflag] = fmincon(cost_func, z0, [], [], [], [], [], [], nonlcon_func, options);


    [X, U, J] = opti_fun(target, rotInertia_arr, mass_out, L/1000);

    % Extract optimal values
    x_star = full([vec(X); vec(U)]);
    f_star = 0;  % No local objective in ATC framework
    
    % Calculate target loss
    % x = reshape(z_opt(1:2*NB*N), 2*NB, N);
    % u = reshape(z_opt(2*NB*N+1:end), NB, N);
    % 
    % Compute performance metrics
    % effort = sum(u(:).^2);
    % smoothness = sum(diff(u,1,2).^2,'all');
    % path_penalty = 0;
    % for k = 1:N-1
    %     p1 = ee_pos(model, x(1:NB,k));
    %     p2 = ee_pos(model, x(1:NB,k+1));
    %     path_penalty = path_penalty + sum((p2-p1).^2);
    % end
    % 
    % r = [effort; smoothness; path_penalty];
    % phi_star = sum((r - target).^2);
    phi_star = full(J);
end

%% Planner subsystem ATC objective function
function f = planner_subsystem_atc_obj(z, target, NB, N, DT, model, alpha)
    % Extract states and inputs from decision vector
    x = reshape(z(1:2*NB*N), 2*NB, N);
    u = reshape(z(2*NB*N+1:end), NB, N);
    
    % Compute performance metrics
    effort = sum(u(:).^2);
    smoothness = sum(diff(u,1,2).^2,'all');
    end_vel = sum(x(NB+1:end,end).^2);
    path_penalty = 0;
    for k = 1:N-1
        p1 = ee_pos(model, x(1:NB,k));
        p2 = ee_pos(model, x(1:NB,k+1));
        path_penalty = path_penalty + sum((p2-p1).^2);
    end
    
    % Original objective function
    J_original = effort + smoothness + end_vel + alpha * path_penalty;
    
    % ATC objective: original objective + target deviation
    r = [effort; smoothness; path_penalty];
    J_target = sum((r - target).^2);
    
    % Final objective is a combination
    f = J_original + 10*J_target;  % Weight the target deviation higher
end

%% Nonlinear constraints for path planning
function [c, ceq] = nonlcon_planner(z, NB, N, DT, model, q0, qf, tau_max)
    x = reshape(z(1:2*NB*N), 2*NB, N);
    u = reshape(z(2*NB*N+1:end), NB, N);
    
    % Dynamic constraints - reduced for computational efficiency
    % Only check every 5 steps to reduce computation
    step_skip = 5;
    n_check_steps = floor((N-1)/step_skip);
    
    ceq_dyn = zeros(2*NB*n_check_steps, 1);
    idx = 0;
    
    for j = 1:n_check_steps
        k = (j-1)*step_skip + 1;
        qk = x(1:NB, k);
        dqk = x(NB+1:end, k);
        ddq = FDab(model, qk, dqk, u(:, k));
        h = x(:, k) + [dqk; ddq]*DT*step_skip - x(:, k+step_skip);
        ceq_dyn(idx+1:idx+2*NB) = h;
        idx = idx + 2*NB;
    end
    
    % Boundary conditions
    ceq_bc = [x(1:NB,1)-q0;
              x(1:NB,end)-qf;
              x(NB+1:end,1);
              x(NB+1:end,end)];
    
    ceq = [ceq_dyn; ceq_bc];
    
    % Input and state bounds - sample for computational efficiency
    sample_indices = 1:5:N;
    c1 = x(1:NB, sample_indices) - pi;
    c2 = -pi - x(1:NB, sample_indices);
    c3 = u(:, sample_indices) - tau_max;
    c4 = -tau_max - u(:, sample_indices);
    c = [c1(:); c2(:); c3(:); c4(:)];
end

%% Forward Dynamics for robot model
function ddq = FDab(model, q, dq, tau)
    NB = length(q);
    
    % Simplified forward dynamics (using inverse dynamics)
    % In a full implementation, use the spatial algebra library for better accuracy
    
    % Approximate mass matrix (simplified)
    M = zeros(NB);
    for i = 1:NB
        M(i,i) = model.I{i}(1,1);  % Diagonal approximation
    end
    
    % Compute Coriolis and gravity terms using ID with zero acceleration
    tau_bias = ID(model, q, dq, zeros(NB,1));
    
    % Solve for accelerations
    ddq = M \ (tau - tau_bias);
end

%% Inverse Dynamics
function tau = ID(model, q, dq, ddq)
    NB = length(q);
    g = [0; 0; -9.81];  % Gravity vector
    
    % Simplified inverse dynamics (placeholder)
    tau = zeros(NB, 1);
    
    % Simplified model: tau = M*ddq + C(q,dq) + G(q)
    for i = 1:NB
        % Approximate terms from the I matrix
        m_i = model.I{i}(4,4);  % Mass from spatial inertia diagonal
        
        % Inertia term
        tau(i) = model.I{i}(1,1) * ddq(i);
        
        % Simple Coriolis approximation
        for j = 1:NB
            if j ~= i
                tau(i) = tau(i) + 0.1 * m_i * dq(j)^2 * sin(q(i)-q(j));
            end
        end
        
        % Simple gravity term
        tau(i) = tau(i) + m_i * 9.81 * cos(q(i));
    end
end

%% Controller Subsystem ATC Optimization
function [x_star, f_star, phi_star, exitflag] = controller_subsystem_atc_opt(target, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max)
    % Inputs:
    %   target - 2×1 vector: [RMSE_target; phase_margin_target]
    %   qd_desired - Desired trajectory (1×N)
    %   I_func, tau_ext_func - Time-varying inertia and disturbance functions
    %   dt, N, q_res, tau_max - Simulation parameters
    % Outputs:
    %   x_star - 3×1 vector of optimal [Kp; Ki; Kd]
    %   f_star - Local objective
    %   phi_star - ATC target loss
    %   exitflag - Optimization exit flag
    
    % Objective function
    objective = @(K) controller_subsystem_atc_obj(K, target, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    
    % Bounds on PID gains
    lb = [0, 0, 0];       % Non-negative gains
    ub = [500, 50, 50];   % Upper bounds to avoid instability
    
    % Initial guess
    initial_guess = [10.0, 1.0, 0.1];
    
    % Optimization options - reduced for multi-run
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'FiniteDifferenceType', 'forward', ...
        'StepTolerance', 1e-3, ...
        'OptimalityTolerance', 1e-3, ...
        'MaxFunctionEvaluations', 100);
    
    % Run optimization
    [K_opt, fval, exitflag] = fmincon(objective, initial_guess, [], [], [], [], lb, ub, @phase_margin_constraint, options);
    
    % Return optimization results
    x_star = K_opt';
    f_star = 0;  % No local objective in ATC framework
    
    % Calculate target loss
    [RMSE, ~, ~] = pid_simulation(K_opt, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    phi_m = phase_margin(K_opt);
    r = [RMSE; phi_m];
    phi_star = sum((r - target).^2);
end

%% Controller subsystem ATC objective function
function f = controller_subsystem_atc_obj(K, target, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max)
    % Simulate controller and compute RMSE
    [RMSE, ~, ~] = pid_simulation(K, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    
    % Compute phase margin
    phi_m = phase_margin(K);
    
    % Target metrics
    r = [RMSE; phi_m];
    
    % ATC objective: target deviation
    f = sum((r - target).^2);
end

%% PID simulation function
function [RMSE, q, error] = pid_simulation(K, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max)
    % Extract PID gains
    Kp = K(1); Ki = K(2); Kd = K(3);
    
    % Initialize state variables
    q = zeros(1, N);
    dq = zeros(1, N);
    integral_error = 0;
    
    % PID control loop
    for k = 3:N
        % Current error
        e = qd_desired(k) - q(k-1);
        
        % Integral term
        integral_error = integral_error + e * dt;
        
        % Derivative term (using backward difference)
        derivative_error = (qd_desired(k) - q(k-1) - (qd_desired(k-1) - q(k-2))) / dt;
        
        % PID control law
        tau = Kp * e + Ki * integral_error + Kd * derivative_error;
        
        % Apply torque limit
        tau = max(min(tau, tau_max), -tau_max);
        
        % System dynamics: I(t) * d²q/dt² = τ - τ_ext
        I_k = I_func(k);
        tau_ext_k = tau_ext_func(k);
        ddq = (tau - tau_ext_k) / I_k;
        
        % Numerical integration (Euler method)
        dq(k) = dq(k-1) + ddq * dt;
        q(k) = q(k-1) + dq(k) * dt;
        
        % Enforce position resolution (quantization)
        q(k) = round(q(k) / q_res) * q_res;
    end
    
    % Compute tracking error and RMSE
    error = qd_desired - q;
    RMSE = sqrt(mean(error.^2));
end

%% Phase margin constraint function
function [c, ceq] = phase_margin_constraint(K)
    phi_m = phase_margin(K);
    
    % Constraint: phi_m <= 45 degrees
    c = phi_m - 45;  % c <= 0
    ceq = [];
end

%% Phase margin calculation
function phi_m = phase_margin(K)
    % Extract PID gains
    Kp = K(1); Ki = K(2); Kd = K(3);
    
    % Use control system toolbox if available, otherwise approximate
    try
        % Define transfer functions
        s = tf('s');
        I_nominal = 1;  % Nominal inertia
        G = 1 / (I_nominal * s^2);
        C = Kp + Ki/s + Kd*s;
        
        % Open-loop transfer function
        L = C * G;
        
        % Compute phase margin
        [~, phi_m, ~] = margin(L);
        
        % Handle invalid results
        if isnan(phi_m) || isinf(phi_m)
            phi_m = 180;
        end
    catch
        % Simple approximation if control system toolbox not available
        if Ki == 0
            phi_m = 90;
        else
            phi_m = atan2(Kp, Ki) * 180/pi;
        end
        
        % Add contribution from derivative term
        phi_m = phi_m + 15 * (Kd / (Kp + 1));
        
        % Limit range
        phi_m = max(min(phi_m, 180), 0);
    end
end

%% Robot Model Creation
function model = create3DoFRobotModel(rotInertia_arr, mass, L)
    % Create 3-DOF robot model with spatial inertia matrices
    NB = 3;
    model.NB = NB;
    model.parent = [0, 1, 2];
    types = {'Rz', 'Rx', 'Rx'};
    
    for i = 1:NB
        model.jtype{i} = types{i};
        model.Xtree{i} = plux(eye(3), [0; 0; L(i)]);
        model.I{i} = mcI(mass(i), [0, 0, L(i)], diag(rotInertia_arr(:,i)));
    end
    
    model.Xtree{NB+1} = plux(eye(3), [0; 0; L(NB)]);
end

%% Spatial Inertia Functions
function I = mcI(m, c, Ic)
    % Spatial inertia for a body with mass m, CoM c, and inertia Ic
    I = zeros(6, 6);
    
    % Mass properties
    I(1:3, 1:3) = Ic;
    I(4:6, 4:6) = m * eye(3);
    
    % Cross-coupling terms for offset CoM
    if norm(c) > 0
        cx = skew(c);
        I(1:3, 4:6) = m * cx;
        I(4:6, 1:3) = m * cx';
        I(1:3, 1:3) = I(1:3, 1:3) + m * cx * cx';
    end
end

function S = skew(v)
    % 3×3 skew-symmetric matrix
    S = [   0, -v(3),  v(2);
          v(3),    0, -v(1);
         -v(2), v(1),    0];
end

function X = plux(R, p)
    % 6×6 spatial transform from rotation R and translation p
    X = [R, zeros(3); skew(p)*R, R];
end

function p = ee_pos(model, q)
    % Compute end-effector position for a given joint configuration
    NB = model.NB;
    
    % Initialize transforms
    X0 = cell(NB, 1);
    
    % Forward kinematics
    for i = 1:NB
        [XJ, ~] = jcalc(model.jtype{i}, q(i));
        if model.parent(i) == 0
            X0{i} = XJ * model.Xtree{i};
        else
            X0{i} = X0{model.parent(i)} * XJ * model.Xtree{i};
        end
    end
    
    % End-effector transform
    X_ee = X0{NB} * model.Xtree{NB+1};
    
    % Extract position
    [~, p] = plux_inv(X_ee);
end

function [XJ, S] = jcalc(type, q)
    % Joint transform and motion subspace
    XJ = eye(6);
    S = zeros(6, 1);
    
    switch type
        case 'Rx'  % Rotation about x-axis
            XJ = [rotx(q), zeros(3); zeros(3), rotx(q)];
            S = [1; 0; 0; 0; 0; 0];
        case 'Ry'  % Rotation about y-axis
            XJ = [roty(q), zeros(3); zeros(3), roty(q)];
            S = [0; 1; 0; 0; 0; 0];
        case 'Rz'  % Rotation about z-axis
            XJ = [rotz(q), zeros(3); zeros(3), rotz(q)];
            S = [0; 0; 1; 0; 0; 0];
    end
end

function R = rotx(theta)
    % Rotation matrix about x-axis
    R = [1, 0, 0;
         0, cos(theta), -sin(theta);
         0, sin(theta), cos(theta)];
end

function R = roty(theta)
    % Rotation matrix about y-axis
    R = [cos(theta), 0, sin(theta);
         0, 1, 0;
         -sin(theta), 0, cos(theta)];
end

function R = rotz(theta)
    % Rotation matrix about z-axis
    R = [cos(theta), -sin(theta), 0;
         sin(theta), cos(theta), 0;
         0, 0, 1];
end

function [R, p] = plux_inv(X)
    % Extract rotation and position from spatial transform
    R = X(1:3, 1:3);
    S = X(4:6, 1:3) * R';
    p = [S(3,2); S(1,3); S(2,1)];
end

%% ATC optimization helper function
function [x_star, f_star, phi_star, exitflag] = atc_opt(atc_obj, target, weight, x0, A, b, Aeq, beq, lb, ub, nonlcon, options, quiet)
    % Define a dummy scalar function of x only 
    function f_atc = atc_obj_dummy(x)  
        [f, phi, phi_sub_star] = atc_obj(x, target);

        % ATC objective function: f + w*target_loss + w*subsys_target_loss*
        f_atc = f + weight(1)*phi + weight(2)*phi_sub_star;
    end

    if quiet
        % Solve with optimization display off
        options.Display = 'off';
        [x_star, f_atc_star, exitflag] = fmincon(@atc_obj_dummy, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    else
        % Solve with default or custom display
        [x_star, f_atc_star, exitflag] = fmincon(@atc_obj_dummy, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    end
  
    % Evaluate again with x* to recover each term of ATC objective function
    [f_star, phi_star, ~] = atc_obj(x_star, target);
end
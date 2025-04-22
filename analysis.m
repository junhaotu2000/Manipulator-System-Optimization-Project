%% Overall Script for Post-Optimality Sensitivity Analysis

clear all; 
close all hidden;
clc;
run setup.m;  % Load required libraries and functions
import casadi.*

%% 1. Mechanical Design Optimization
[rotInertia_arr, mass_out, L] = MechOptimization();
fprintf('Mechanical Optimization Complete.\n');
disp('Rotational Inertia Array (each column corresponds to [I_x; I_y; I_z]):');
disp(rotInertia_arr);
disp('Mass Vector (link masses in kg):');
disp(mass_out);

%% 2. Nominal Robot Planning Optimization via CasADi
% Create the optimization function
opti_fun = make_arm_opti();

% --- Capture the outputs as they are returned by your optimizer ---
% For example, if your function returns two outputs:
[x_nom, u_nom] = opti_fun(rotInertia_arr', mass_out, L/1000);

% If your optimizer also returns dual variables, use:
% [x_nom, u_nom, lambda_nom] = opti_fun(rotInertia_arr', mass_out, L/1000);

% Convert solutions from CasADi DM to full MATLAB arrays
X_opt_nom = full(x_nom);
U_opt_nom = full(u_nom);

% Display the nominal results
fprintf('Nominal Robot Planning Optimization Complete.\n');

% If dual variables were returned, display them:
if exist('lambda_nom', 'var')
    fprintf('Dual Variables (Lagrange Multipliers):\n');
    disp(full(lambda_nom));
else
    fprintf('No dual variable information available from the solver.\n');
end

% Compute the nominal objective value
% NOTE: Ensure you have implemented computeObjective to accept the outputs from opti_fun.
J_nom = computeObjective(x_nom, u_nom);
fprintf('Nominal Objective Value: %f\n', J_nom);

%% 3. Post-Optimality Sensitivity Analysis using Finite Difference

% 3.1 Sensitivity with Respect to Link Length (L)
epsilon = 0.01;  % 1% perturbation
L_perturbed = L * (1 + epsilon);

% Re-run the optimization with perturbed L
[x_pert_L, u_pert_L] = opti_fun(rotInertia_arr', mass_out, L_perturbed/1000);
J_pert_L = computeObjective(x_pert_L, u_pert_L);

% Finite-difference sensitivity approximation:
sensitivity_L = (J_pert_L - J_nom) / (epsilon * L);
fprintf('Sensitivity of the cost with respect to link length L: %f\n', sensitivity_L);

% 3.2 Sensitivity with Respect to Link Masses (mass_out)
epsilon_mass = 0.01;  % 1% perturbation in mass
mass_perturbed = mass_out * (1 + epsilon_mass);

% Re-run the optimization with perturbed mass
[x_pert_mass, u_pert_mass] = opti_fun(rotInertia_arr', mass_perturbed, L/1000);
J_pert_mass = computeObjective(x_pert_mass, u_pert_mass);

% Calculate sensitivity for each mass component (assuming elementwise effect)
sensitivity_mass = (J_pert_mass - J_nom) ./ (epsilon_mass .* mass_out);
fprintf('Sensitivity of the cost with respect to mass components:\n');
disp(sensitivity_mass);

% (Optionally, you can add similar perturbations for the inertia parameters.)
%% 5. (Optional) Plotting the Nominal Trajectories
time_grid = linspace(0, 2, size(X_opt_nom,2));  % Assuming a total duration of 2 seconds

figure('Name','Nominal Robot States & Inputs','NumberTitle','off');
subplot(1,3,1); % Joint Angles
plot(time_grid, X_opt_nom(1,:), 'LineWidth', 1.5); hold on;
plot(time_grid, X_opt_nom(2,:), 'LineWidth', 1.5);
plot(time_grid, X_opt_nom(3,:), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
title('Nominal Joint Angles'); grid on;

subplot(1,3,2); % (If X_opt_nom includes velocities; adjust as needed)
% For example:
% plot(time_grid, X_opt_nom(4,:), 'LineWidth', 1.5); hold on;
% plot(time_grid, X_opt_nom(5,:), 'LineWidth', 1.5);
% plot(time_grid, X_opt_nom(6,:), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (rad/s)');
title('Nominal Joint Velocities'); grid on;

subplot(1,3,3); % Joint Torques
plot(time_grid, U_opt_nom(1,:), 'LineWidth', 1.5); hold on;
plot(time_grid, U_opt_nom(2,:), 'LineWidth', 1.5);
plot(time_grid, U_opt_nom(3,:), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Torque (Nm)');
title('Nominal Joint Torques'); grid on;

%% 4. Controller Optimization
% The controller optimization is run for each joint separately.
% Here, we assume that run_pid_optimization uses the optimal control (torque) and joint angle trajectories.
I_const = 0.02 * ones(size(U_opt_nom,2), 1);  % Example fixed inertia vector

% Joint 1
[RMSE_Final_1, q1] = run_pid_optimization(I_const, U_opt_nom(1,:), X_opt_nom(1,:));

% Joint 2
[RMSE_Final_2, q2] = run_pid_optimization(I_const, U_opt_nom(2,:), X_opt_nom(2,:));

% Joint 3
[RMSE_Final_3, q3] = run_pid_optimization(I_const, U_opt_nom(3,:), X_opt_nom(3,:));

% Combine controller results for reporting
RMSE_Final = [RMSE_Final_1; RMSE_Final_2; RMSE_Final_3];
q = [q1; q2; q3];
fprintf('Controller Optimization Complete.\n');


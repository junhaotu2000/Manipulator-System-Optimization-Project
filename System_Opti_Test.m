%% 
% Enhanced Manipulator System Optimization using ATC (Analytical Target Cascading)
% This script integrates three subsystems with ATC coordination and runs multiple optimizations
% 1. Mechanical subsystem: GA-based link optimization with structural constraints
% 2. Path planning subsystem: Trajectory optimization with dynamics and path constraints
% 3. Control subsystem: PID controller tuning with time-varying dynamics
%

clear all;
close all;
clc;

% Global settings
visualize_results = false;  % Set to true to enable visualization
num_optimization_runs = 5;  % Number of times to run the full optimization

%% Top-level system targets
% System-level targets for coupling variables
target_RMSE = 0.05;            % Tracking error target
target_effort = 50;            % Effort minimization target (reduced)
target_smoothness = 20;        % Trajectory smoothness target (increased)
target_path_penalty = 0.5;     % Path deviation penalty target (reduced)
target_phase_margin = 40;      % Phase margin target (degrees) (increased for stability)

% Combined target vector
system_targets = [target_RMSE; target_effort; target_smoothness; target_path_penalty; target_phase_margin];

% Display system targets
fprintf('\n============= System-level Targets =============\n');
fprintf('Tracking Error (RMSE):      %.4f\n', target_RMSE);
fprintf('Control Effort:             %.4f\n', target_effort);
fprintf('Trajectory Smoothness:      %.4f\n', target_smoothness);
fprintf('Path Penalty:               %.4f\n', target_path_penalty);
fprintf('Phase Margin:               %.4f degrees\n', target_phase_margin);
fprintf('==============================================\n\n');

%% Run multiple ATC optimizations and store results
% Initialize storage for results
all_x_star = zeros(11, num_optimization_runs);
all_f_star = zeros(1, num_optimization_runs);
all_phi_star = zeros(1, num_optimization_runs);
all_exitflags = zeros(1, num_optimization_runs);
all_mech_params = zeros(9, num_optimization_runs);  % D, t, L
all_control_params = zeros(3, num_optimization_runs);  % Kp, Ki, Kd

for run_idx = 1:num_optimization_runs
    fprintf('\n============= Starting Optimization Run %d/%d =============\n', run_idx, num_optimization_runs);
    
    % Run system-level ATC optimization
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
    
    % Run mechanical subsystem optimization
    fprintf('Solving mechanical subsystem with GA optimization...\n');
    mech_x_star = mechanical_subsystem_atc_opt(mech_targets);
    all_mech_params(:, run_idx) = mech_x_star;
    
    % Extract link properties
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
    
    % Run path planning optimization (with smoothness improvements)
    fprintf('Solving path planning subsystem with improved dynamics constraints...\n');
    plan_x_star = planner_subsystem_atc_opt(plan_targets, rotInertia_arr, mass_out, L);
    
    % Extract trajectory for controller
    x = reshape(plan_x_star(1:6*N), 6, N);
    q_planned = x(1:3,:);
    
    % Apply smoothing to the trajectory for better control stability
    q_smoothed = q_planned;
    for i = 1:3
        q_smoothed(i,:) = smoothdata(q_planned(i,:), 'gaussian', 5);
    end
    
    % Define time-varying dynamics for controller
    t_vec = linspace(0, T, N);
    I_func = @(k) 1 + 0.1*sin(t_vec(min(k, N)));  % Reduced variation for stability
    tau_ext_func = @(k) 0.05*t_vec(min(k, N));    % Reduced disturbance
    qd_desired = q_smoothed(1,:);                 % Use smoothed trajectory
    
    % Run controller optimization with stability improvements
    fprintf('Solving controller subsystem with improved stability constraints...\n');
    ctrl_x_star = controller_subsystem_atc_opt(ctrl_targets, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    all_control_params(:, run_idx) = ctrl_x_star;
    
    % Display results for this run
    fprintf('----------- Run %d Results -----------\n', run_idx);
    fprintf('System objective (f_star): %.4f\n', f_star);
    fprintf('Target deviation (phi_star): %.4f\n', phi_star);
    fprintf('Exit flag: %d\n', exitflag);
    fprintf('Link diameters (mm): [%.2f, %.2f, %.2f]\n', D(1), D(2), D(3));
    fprintf('Link thicknesses (mm): [%.2f, %.2f, %.2f]\n', t(1), t(2), t(3));
    fprintf('Link lengths (mm): [%.2f, %.2f, %.2f]\n', L(1), L(2), L(3));
    fprintf('Total mass (kg): %.4f\n', sum(mass_out));
    fprintf('Controller gains - Kp: %.4f, Ki: %.4f, Kd: %.4f\n', ctrl_x_star(1), ctrl_x_star(2), ctrl_x_star(3));
    fprintf('--------------------------------------\n');
    
    % Optional visualization for the final run
    if visualize_results && run_idx == num_optimization_runs
        visualize_system_solution(t_vec, q_smoothed, x, plan_x_star, rotInertia_arr, mass_out, L, ctrl_x_star, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);
    end
end

%% Tabulate and display all optimization results
fprintf('\n============= Summary of %d Optimization Runs =============\n', num_optimization_runs);

% System-level objectives table
fprintf('\nSystem-Level Objectives:\n');
results_table = table((1:num_optimization_runs)', all_f_star', all_phi_star', all_exitflags', ...
                      'VariableNames', {'Run', 'ObjectiveValue', 'TargetDeviation', 'ExitFlag'});
disp(results_table);

% System-level design variables table
fprintf('\nSystem-Level Design Variables (x_star):\n');
var_names = {'Run', 'I_trans_1', 'I_trans_2', 'I_trans_3', 'Mass_1', 'Mass_2', 'Mass_3', ...
             'Effort', 'Smoothness', 'PathPenalty', 'RMSE', 'PhaseMargin'};
design_table = array2table([(1:num_optimization_runs)', all_x_star'], 'VariableNames', var_names);
disp(design_table);

% Mechanical subsystem results
fprintf('\nMechanical Subsystem Results:\n');
mech_table = table((1:num_optimization_runs)', ...
                  all_mech_params(1,:)', all_mech_params(2,:)', all_mech_params(3,:)', ...
                  all_mech_params(4,:)', all_mech_params(5,:)', all_mech_params(6,:)', ...
                  all_mech_params(7,:)', all_mech_params(8,:)', all_mech_params(9,:)', ...
                  'VariableNames', {'Run', 'D1', 'D2', 'D3', 't1', 't2', 't3', 'L1', 'L2', 'L3'});
disp(mech_table);

% Controller subsystem results
fprintf('\nController Subsystem Results:\n');
ctrl_table = table((1:num_optimization_runs)', all_control_params(1,:)', all_control_params(2,:)', all_control_params(3,:)', ...
                   'VariableNames', {'Run', 'Kp', 'Ki', 'Kd'});
disp(ctrl_table);

% Calculate statistics across runs
fprintf('\nStatistics Across %d Runs:\n', num_optimization_runs);
fprintf('Mean system objective: %.4f ± %.4f\n', mean(all_f_star), std(all_f_star));
fprintf('Mean target deviation: %.4f ± %.4f\n', mean(all_phi_star), std(all_phi_star));
fprintf('Successful runs: %d/%d\n', sum(all_exitflags > 0), num_optimization_runs);

% Find the best run (minimum f_star)
[~, best_idx] = min(all_f_star);
fprintf('\nBest Design (Run %d):\n', best_idx);
fprintf('System objective (f_star): %.4f\n', all_f_star(best_idx));
fprintf('Target deviation (phi_star): %.4f\n', all_phi_star(best_idx));
fprintf('Link diameters (mm): [%.2f, %.2f, %.2f]\n', all_mech_params(1:3, best_idx)');
fprintf('Link thicknesses (mm): [%.2f, %.2f, %.2f]\n', all_mech_params(4:6, best_idx)');
fprintf('Link lengths (mm): [%.2f, %.2f, %.2f]\n', all_mech_params(7:9, best_idx)');
fprintf('Controller gains - Kp: %.4f, Ki: %.4f, Kd: %.4f\n', all_control_params(:, best_idx)');
fprintf('================================================================\n');

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
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                          'MaxFunctionEvaluations', 10000, 'MaxIterations', 100);
    
    % Optimize using ATC helper function
    [x_star, f_star, phi_star, exitflag] = atc_opt(@manipulator_system_atc_obj, target, weight, x0, A, b, Aeq, beq, lb, ub, nonlcon, options, false);
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
    if all(target < 1e-10)
        % Fast evaluation for system level
        x_star = (lb + ub)' / 2;  % Middle point for default solution
        [f_star, phi_star, ~] = mechanical_subsystem_atc_obj(x_star, target);
        exitflag = 1;
        return;
    end
    
    % Build initial population (similar to MechOptimization.m)
    popSize = 100;  % Reduced population for faster convergence in ATC
    nCrit = 20;     % Number of individuals at midpoint
    
    initialPop = zeros(popSize, nCont);
    for i = 1:nCrit
        initialPop(i,:) = (lb + ub) / 2;
    end
    for i = (nCrit+1):popSize
        initialPop(i,:) = lb + rand(1, nCont) .* (ub - lb);
    end
    
    % GA options
    optionsGA = optimoptions('ga', ...
        'Display', 'off', ...           % Turn off display for ATC integration
        'PopulationSize', popSize, ...
        'MaxGenerations', 50, ...       % Reduced generations for ATC
        'InitialPopulationMatrix', initialPop, ...
        'UseParallel', false);          % Disable parallel for integration
    
    % Run GA optimization
    [x_opt, fval_opt, exitflag] = ga(atc_obj, nCont, [], [], [], [], lb, ub, @nonlcon_mech, [], optionsGA);
    
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
function [x_star, f_star, phi_star, exitflag] = planner_subsystem_atc_opt(target, rotInertia_arr, mass_out, L)
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
    
    % Create robot model
    model = create3DoFRobotModel(rotInertia_arr, mass_out, L/1000);
    
    % Parameters
    NB = model.NB;
    N = 100;
    T = 2.0;
    DT = T/N;
    q_init = [0; 0; 0];
    q_goal = [pi/2; pi/4; 0];
    tau_max = 15;
    alpha = 5;
    
    % Initial guess
    z0 = zeros((2*NB + NB)*N, 1);
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 10000, ...
        'OptimalityTolerance', 1e-4, ...
        'StepTolerance', 1e-4);
    
    % Create objective function with embedded target
    cost_func = @(z) planner_subsystem_atc_obj(z, target, NB, N, DT, model, alpha);
    
    % Nonlinear constraints
    nonlcon_func = @(z) nonlcon_planner(z, NB, N, DT, model, q_init, q_goal, tau_max);
    
    % Run optimization
    [z_opt, fval, exitflag] = fmincon(cost_func, z0, [], [], [], [], [], [], nonlcon_func, options);
    
    % Extract optimal values
    x_star = z_opt;
    f_star = 0;  % No local objective in ATC framework
    
    % Calculate target loss
    x = reshape(z_opt(1:2*NB*N), 2*NB, N);
    u = reshape(z_opt(2*NB*N+1:end), NB, N);
    
    % Compute performance metrics
    effort = sum(u(:).^2);
    smoothness = sum(diff(u,1,2).^2,'all');
    path_penalty = 0;
    for k = 1:N-1
        p1 = ee_pos(model, x(1:NB,k));
        p2 = ee_pos(model, x(1:NB,k+1));
        path_penalty = path_penalty + sum((p2-p1).^2);
    end
    
    r = [effort; smoothness; path_penalty];
    phi_star = sum((r - target).^2);
end

%% Planner subsystem ATC objective function
function f = planner_subsystem_atc_obj(z, target, NB, N, DT, model, alpha)
    % Extract states and inputs from decision vector
    x = reshape(z(1:2*NB*N), 2*NB, N);
    u = reshape(z(2*NB*N+1:end), NB, N);
    
    % Compute performance metrics with enhanced smoothness penalties
    effort = sum(u(:).^2);
    
    % Enhanced smoothness metric - penalize both torque and state derivatives
    smoothness_torque = sum(diff(u,1,2).^2,'all');
    smoothness_vel = sum(diff(x(NB+1:end,:),1,2).^2,'all');
    smoothness_pos = sum(diff(diff(x(1:NB,:),1,2),1,2).^2,'all');
    smoothness = smoothness_torque + 0.5*smoothness_vel + 2.0*smoothness_pos;
    
    % End velocity constraint (stop at goal)
    end_vel = sum(x(NB+1:end,end).^2);
    
    % Path penalty - penalize jerky movements
    path_penalty = 0;
    for k = 1:N-2
        p1 = ee_pos(model, x(1:NB,k));
        p2 = ee_pos(model, x(1:NB,k+1));
        p3 = ee_pos(model, x(1:NB,k+2));
        
        % Penalize both distance and change in velocity (jerk)
        path_diff1 = p2 - p1;
        path_diff2 = p3 - p2;
        path_penalty = path_penalty + sum((path_diff2 - path_diff1).^2);
    end
    
    % Enhanced objective function with stronger regularization
    J_original = effort + 5.0*smoothness + 10.0*end_vel + alpha * path_penalty;
    
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
    
    % Dynamic constraints
    ceq_dyn = zeros(2*NB*(N-1),1);
    idx = 0;
    for k = 1:N-1
        qk = x(1:NB,k);
        dqk = x(NB+1:end,k);
        ddq = FDab(model, qk, dqk, u(:,k));
        h = x(:,k) + [dqk; ddq]*DT - x(:,k+1);
        ceq_dyn(idx+1:idx+2*NB) = h;
        idx = idx + 2*NB;
    end
    
    % Boundary conditions
    ceq_bc = [x(1:NB,1)-q0;
              x(1:NB,end)-qf;
              x(NB+1:end,1);
              x(NB+1:end,end)];
    
    ceq = [ceq_dyn; ceq_bc];
    
    % Input and state bounds
    c = [x(1:NB,:) - pi;
         -pi - x(1:NB,:);
         u - tau_max;
         -tau_max - u];
    c = c(:);
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
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'FiniteDifferenceType', 'central', ...
        'StepTolerance', 1e-4, ...
        'OptimalityTolerance', 1e-4, ...
        'MaxFunctionEvaluations', 300);
    
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
        % Define transfer functions with improved plant model
        s = tf('s');
        I_nominal = 1;  % Nominal inertia
        
        % More accurate plant model with damping
        zeta = 0.05;  % Small damping factor for stability
        omega_n = 10; % Natural frequency
        G = 1 / (I_nominal * s^2 + 2*zeta*omega_n*I_nominal*s);
        
        % Standard PID controller
        C = Kp + Ki/s + Kd*s/(1 + s/100);  % Filter on derivative term
        
        % Open-loop transfer function
        L = C * G;
        
        % Compute phase margin with warning suppression
        warning('off', 'Control:analysis:MarginUnstable');
        [~, phi_m, ~] = margin(L);
        warning('on', 'Control:analysis:MarginUnstable');
        
        % Handle invalid results
        if isnan(phi_m) || isinf(phi_m)
            phi_m = 0;  % Assume poor phase margin for invalid results
        end
    catch
        % Improved approximation if control system toolbox not available
        % Based on standard formulas for second-order systems with PID
        
        % Simplified approximation
        if Ki == 0
            if Kd == 0
                phi_m = 90;  % P controller
            else
                phi_m = 90 + atan2(Kd, Kp) * 180/pi;  % PD controller
            end
        else
            if Kd == 0
                phi_m = atan2(Kp, Ki) * 180/pi;  % PI controller
            else
                % Full PID
                % Calculate rough phase margin based on gain and time constants
                Ti = Kp/Ki;  % Integral time constant
                Td = Kd/Kp;  % Derivative time constant
                
                % Approximate phase margin based on time constants
                if Ti > 4*Td
                    phi_m = 45 + 25*log10(Td) + 25*log10(Ti);
                else
                    phi_m = 30 + 30*log10(Td);
                end
            end
        end
        
        % Limit range to realistic values for control systems
        phi_m = max(min(phi_m, 90), 0);
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
clc; clear; close all;

%% Simulation Parameters
T = 10;          % Total simulation time (s)
dt = 1e-3;       % Temporal resolution (1 ms)
N = T/dt;        % Number of discrete time steps
q0 = 0;          % Initial position (rad)
dq0 = 0;         % Initial velocity (rad/s)
q_res = 0.18;    % Angular position resolution (rad)
tau_max = 300;   % Actuator torque limit (Nm)

% Time vector
t = linspace(0, T, N);

% Desired trajectory (example: sinusoidal)
%qd_desired = sin(0.5 * t);
qd_desired = cos(t) + 0.25*sin(3*t);



%% Time-Varying Inertia and External Torque
I_func = 1 + 0.5*sin(t);  % Example: Fluctuating inertia
%I_func = -0.3;



tau_ext_func = t; % Example: External torque variations
%tau_ext_func = 0.7;



%% PID Controller Simulation Function
function [RMSE, q, error] = pid_simulation(K, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max)
    Kp = K(1); Ki = K(2); Kd = K(3);
    
    % Initialize state variables
    q = zeros(1, N);
    dq = zeros(1, N);
    integral_error = 0;
    
    % PID control loop with fixed indexing
    for k = 3:N
        error = qd_desired(k) - q(k-1);
        integral_error = integral_error + error * dt;
        
        % Compute derivative of error correctly
        derivative_error = (qd_desired(k) - q(k-1) - (qd_desired(k-1) - q(k-2))) / dt;
        
        % PID control law
        tau = Kp * error + Ki * integral_error + Kd * derivative_error;
        
        % Apply torque limit
        tau = max(min(tau, tau_max), -tau_max);
        
        % System dynamics: I(t) * d²q/dt² = τ - τ_ext
        I_k = I_func(k);
        tau_ext_k = tau_ext_func(k);
        ddq = (tau - tau_ext_k) / I_k;
        
        % Numerical integration (Euler method)
        dq(k) = dq(k-1) + ddq * dt;
        q(k) = q(k-1) + dq(k) * dt;
        
        % Enforce angular position resolution constraint
        %q(k) = round(q(k) / q_res) * q_res;
    end
    
    % Compute RMSE
    error = qd_desired - q;
    RMSE = sqrt((1/N) * sum(error.^2))
end

%% Optimization Setup
% Define objective function
objective = @(K) pid_simulation(K, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);

% Bounds on PID gains (Only non-negative)
lb = [0, 0, 0]; % Kp, Ki, Kd must be >= 0
ub = [500, 50, 50]; % Reasonable upper bounds to avoid instability

% Initial Guess (try different values)
initial_guess = [5.0, 0.5, 0.05];

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'FiniteDifferenceType', 'central', ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'MaxFunctionEvaluations', 1000);

% Run optimization
K_opt = fmincon(objective, initial_guess, [], [], [], [], lb, ub, @phase_margin_constraint, options);

% Extract optimized gains
Kp_opt = K_opt(1);
Ki_opt = K_opt(2);
Kd_opt = K_opt(3);

fprintf('Optimized PID Gains: Kp=%.3f, Ki=%.3f, Kd=%.3f\n', Kp_opt, Ki_opt, Kd_opt);

% Run simulation with optimized gains
[RMSE_final, q_final, error_final] = pid_simulation(K_opt, qd_desired, I_func, tau_ext_func, dt, N, q_res, tau_max);

%% Plotting
figure;
subplot(2,1,1);
plot(t, qd_desired, 'r--', 'LineWidth', 1.5); hold on;
plot(t, q_final, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Position (rad)');
legend('Desired Position', 'Actual Position');
title('Time Response of q(t)'); grid on;

subplot(2,1,2);
plot(t, error_final, 'k', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Tracking Error (rad)');
title('Tracking Error e(t) = q_d - q'); grid on;

%% Phase Margin Constraint Function
function [c, ceq] = phase_margin_constraint(K)
    Kp = K(1); Ki = K(2); Kd = K(3);
    
    % Define transfer function G(s) = 1 / (I * s^2)
    s = tf('s');
    I_nominal = 1; % Assuming nominal inertia for phase margin analysis
    G = 1 / (I_nominal * s^2);
    
    % PID Controller Transfer Function
    C = Kp + Ki/s + Kd * s;
    
    % Open-loop transfer function
    L = C * G;
    
    % Compute phase margin
    [~, phi_m, ~] = margin(L);
    
    % Ensure valid values
    if isnan(phi_m) || isinf(phi_m)
        phi_m = 180; % Assign large value to indicate instability
    end
    
    % Constraint: phi_m <= 45 degrees
    c = phi_m - 45;  % Constraint should be <= 0
    ceq = [];
    
    % Debugging: Display phase margin
    fprintf('Phase Margin: %.2f° (Constraint: <= 45°)\n', phi_m);
end

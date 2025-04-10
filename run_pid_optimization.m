function [RMSE_final, q] = run_pid_optimization(I_func, tau_ext_func, qd_desired)

    %% Simulation Parameters
    T = 2;           % Total simulation time (s)
    dt = 2/100;        % Temporal resolution (1 ms)
    N = T/dt;         % Number of discrete time steps
    q_res = 0.18;     % Angular position resolution (rad)
    tau_max = 300;    % Actuator torque limit (Nm)

    if nargin < 3 || isempty(qd_desired)
        qd_desired = cos(t) + 0.25*sin(3*t);  % Default desired trajectory
    end
    if nargin < 2 || isempty(tau_ext_func)
        tau_ext_func = t;  % Default external torque
    end
    if nargin < 1 || isempty(I_func)
        I_func = 1 + 0.5*sin(t);  % Default time-varying inertia
    end

    %% PID Controller Simulation Function
    function [RMSE, q, error] = pid_simulation(K)
        Kp = K(1); Ki = K(2); Kd = K(3);
        
        % Initialize state variables
        q = zeros(1, N);
        dq = zeros(1, N);
        integral_error = 0;
        
        for k = 3:N
            error = qd_desired(k) - q(k-1);
            integral_error = integral_error + error * dt;
            derivative_error = (qd_desired(k) - q(k-1) - (qd_desired(k-1) - q(k-2))) / dt;

            tau = Kp * error + Ki * integral_error + Kd * derivative_error;
            tau = max(min(tau, tau_max), -tau_max);

            I_k = I_func(k);
            tau_ext_k = tau_ext_func(k);
            ddq = (tau - tau_ext_k) / I_k;

            dq(k) = dq(k-1) + ddq * dt;
            q(k) = q(k-1) + dq(k) * dt;
        end

        error = qd_desired - q;
        RMSE = sqrt((1/N) * sum(error.^2));
    end

    %% Phase Margin Constraint Function
    function [c, ceq] = phase_margin_constraint(K)
        Kp = K(1); Ki = K(2); Kd = K(3);
        s = tf('s');
        I_nominal = 1;
        G = 1 / (I_nominal * s^2);
        C = Kp + Ki/s + Kd * s;
        L = C * G;
        [~, phi_m, ~] = margin(L);
        if isnan(phi_m) || isinf(phi_m)
            phi_m = 180;
        end
        c = phi_m - 45;
        ceq = [];
    end

    %% Optimization Setup
    objective = @(K) pid_simulation(K);
    lb = [0, 0, 0];
    ub = [500, 50, 50];
    initial_guess = [5.0, 0.5, 0.05];

    options = optimoptions('fmincon', ...
        'Display', 'none', ...
        'Algorithm', 'sqp', ...
        'FiniteDifferenceType', 'central', ...
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'MaxFunctionEvaluations', 1000);

    % Run optimization
    K_opt = fmincon(objective, initial_guess, [], [], [], [], lb, ub, @phase_margin_constraint, options);

    % Final simulation with optimal gains
    [RMSE_final, ~, ~] = pid_simulation(K_opt);
end

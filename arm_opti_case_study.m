function analyze_constraints()
    clear; clc; close all;
    run setup.m
    import casadi.*
    
    % 1) Robot model and simulation parameters
    model = autoTree(3, 1, pi/2, 1); 
    % Add an additional end-effector offset
    model.Xtree{end+1} = plux( eye(3), [0.1, 0, 0] );
    NB = model.NB; 
    disp(['Robot has NB = ', num2str(NB), ' joints.']);

    T = 2.0;     % Total duration (2 s)
    N = 100;     % Number of discretization steps
    DT = T / N;
    time_grid = linspace(0, T, N);

    % Define multiple sets of constraints to test, stored in a struct array
    test_cases = [
        struct('name','Case1 - AllConstraints', ...
               'max_torque_on',true, ...
               'joint_limit_on',true, ...
               'final_vel_zero',true), ...
        struct('name','Case2 - NoTorqueLimit', ...
               'max_torque_on',false, ...
               'joint_limit_on',true, ...
               'final_vel_zero',true), ...
        struct('name','Case3 - NoJointLimit', ...
               'max_torque_on',true, ...
               'joint_limit_on',false, ...
               'final_vel_zero',true), ...
        struct('name','Case4 - FreeFinalVelocity', ...
               'max_torque_on',true, ...
               'joint_limit_on',true, ...
               'final_vel_zero',false)
    ];

    % Store solutions for subsequent plotting and comparison
    solutions = cell(numel(test_cases),1);

    % 2) Iterate over each constraint combination and solve the optimal control problem
    for iCase = 1:numel(test_cases)
        disp(['=== Solving ', test_cases(iCase).name, ' ===']);
        
        % Create a new Opti object
        opti = casadi.Opti();

        % Decision variables
        x = opti.variable(NB*2, N);  % [q; dq]
        u = opti.variable(NB,   N);  % torque

        % Dynamics constraints (Explicit Euler)
        for k = 1:(N-1)
            xk  = x(:,k);
            uk  = u(:,k);

            qk  = xk(1:NB);
            dqk = xk(NB+1:end);

            ddqk = FDab(model, qk, dqk, uk);
            x_next = xk + [dqk; ddqk]*DT;
            opti.subject_to( x(:,k+1) == x_next );
        end

        % Boundary conditions
        q_init = [0;  0;  0]; 
        q_goal = [pi/2;  0;  pi/2];
        opti.subject_to( x(1:NB,1)   == q_init );
        opti.subject_to( x(1:NB,end) == q_goal );
        
        % Optional: whether initial & final velocity is zero
        if test_cases(iCase).final_vel_zero
            opti.subject_to( x(NB+1:end,1)   == 0 );
            opti.subject_to( x(NB+1:end,end) == 0 );
        else
            % Only enforce zero initial velocity
            opti.subject_to( x(NB+1:end,1)   == 0 );
            % No constraint on final velocity
        end

        % Physical constraints
        if test_cases(iCase).max_torque_on
            max_torque = 15; 
            opti.subject_to( -max_torque <= u <= max_torque );
        end

        if test_cases(iCase).joint_limit_on
            opti.subject_to( -pi <= x(1:NB,:) <= pi );
        end

        % Objective function (multi-objective + end-effector path length)
        control_effort = sum(sum( u.^2 ));
        smoothness_penalty = sum(sum( diff(u,1,2).^2 ));
        end_vel_penalty = sum( x(NB+1:end, end).^2 );
        
        % End-effector path penalty
        alpha_path = 5;
        end_effector_path_length = 0;
        for k = 1 : (N-1)
            p_k   = end_effector_position(model, x(1:NB, k));
            p_kp1 = end_effector_position(model, x(1:NB, k+1));
            end_effector_path_length = end_effector_path_length + sum1( (p_kp1 - p_k).^2 );
        end

        % Combine all cost terms
        J = control_effort + smoothness_penalty + alpha_path * end_effector_path_length;
        
        % If we still want smaller final velocity, we can add an end_vel_penalty
        % Here is an option: only add extra penalty when final_vel_zero=false
        if ~test_cases(iCase).final_vel_zero
            J = J + 10 * end_vel_penalty;
        end
        
        opti.minimize(J);

        % Set the solver
        opti.solver('ipopt', struct('ipopt', struct('print_level', 0)));

        % Initial guess
        opti.set_initial(x, 0);
        opti.set_initial(u, 0);

        % Solve
        sol = opti.solve();

        % Extract the results
        X_opt = sol.value(x);
        U_opt = sol.value(u);

        solutions{iCase}.name = test_cases(iCase).name;
        solutions{iCase}.Q     = X_opt(1:NB,:);
        solutions{iCase}.dQ    = X_opt(NB+1:end,:);
        solutions{iCase}.U     = U_opt;
        solutions{iCase}.cost  = sol.value(J);
    end

    % 3) Results comparison and visualization
    plot_comparison(time_grid, solutions);
end

%% =============== Utility function: end-effector position calculation ===============
function p_ee = end_effector_position(model, q_sym)
    import casadi.*
    NB_ = model.NB;
    Xup = cell(NB_,1);
    X0  = cell(NB_,1);  
    parent = model.parent;
    
    for i = 1 : NB_
        [XJ,~] = jcalc(model.jtype{i}, q_sym(i));
        Xup{i} = XJ * model.Xtree{i};
        if parent(i) == 0
            X0{i} = Xup{i};
        else
            X0{i} = X0{ parent(i) } * Xup{i};
        end
    end
    
    % Coordinate of the end-effector (NB_+1)
    X_ee = X0{NB_} * model.Xtree{NB_+1};
    [~, p_vec] = plux_inv(X_ee);
    p_ee = p_vec;  % 3x1
end

function [R, p] = plux_inv(X)
    R = X(1:3, 1:3);
    S = X(4:6, 1:3) * R'; 
    p = [ S(3,2); S(1,3); S(2,1) ]; 
end

%% =============== Utility function: plot comparison ===============
function plot_comparison(time_grid, solutions)
    % Assume the number of joints is the same for all solutions
    NB = size(solutions{1}.Q, 1);
    
    % --- First row: joint angles ---
    for i = 1:NB
        subplot(3,3, i);
        hold on; grid on;
        plot(time_grid, solutions{1}.Q(i,:), 'LineWidth',1.5, 'DisplayName','Case1');
        plot(time_grid, solutions{2}.Q(i,:), 'LineWidth',1.5, 'DisplayName','Case2');
        plot(time_grid, solutions{3}.Q(i,:), 'LineWidth',1.5, 'DisplayName','Case3');
        plot(time_grid, solutions{4}.Q(i,:), 'LineWidth',1.5, 'DisplayName','Case4');
        title(['Angle - Joint ', num2str(i)]);
        xlabel('Time (s)');
        ylabel('Angle (rad)');
        if i == 1
            legend('Location','northoutside','Orientation','horizontal');
        end
    end
    
    % --- Second row: joint angular velocities ---
    for i = 1:NB
        subplot(3,3, NB + i);
        hold on; grid on;
        plot(time_grid, solutions{1}.dQ(i,:), 'LineWidth',1.5, 'DisplayName','Case1');
        plot(time_grid, solutions{2}.dQ(i,:), 'LineWidth',1.5, 'DisplayName','Case2');
        plot(time_grid, solutions{3}.dQ(i,:), 'LineWidth',1.5, 'DisplayName','Case3');
        plot(time_grid, solutions{4}.dQ(i,:), 'LineWidth',1.5, 'DisplayName','Case4');
        title(['Angular Velocity - Joint ', num2str(i)]);
        xlabel('Time (s)');
        ylabel('Velocity (rad/s)');
    end
    
    % --- Third row: torques ---
    for i = 1:NB
        subplot(3,3, 2*NB + i);
        hold on; grid on;
        plot(time_grid, solutions{1}.U(i,:), 'LineWidth',1.5, 'DisplayName','Case1');
        plot(time_grid, solutions{2}.U(i,:), 'LineWidth',1.5, 'DisplayName','Case2');
        plot(time_grid, solutions{3}.U(i,:), 'LineWidth',1.5, 'DisplayName','Case3');
        plot(time_grid, solutions{4}.U(i,:), 'LineWidth',1.5, 'DisplayName','Case4');
        title(['Torque - Joint ', num2str(i)]);
        xlabel('Time (s)');
        ylabel('Torque (Nm)');
    end
    
    % Display the final cost of each case
    disp('--- Final Costs for Each Case ---');
    for i = 1:numel(solutions)
        disp([solutions{i}.name, '   Cost = ', num2str(solutions{i}.cost)]);
    end
end

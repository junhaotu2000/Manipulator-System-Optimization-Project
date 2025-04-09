function [time_grid, I_eff, torque_traj, angle_traj] = PlannerFunction(rotInertia, mass)
    % PlannerFunction Optimizes the motion of a 3-DOF robot arm and outputs:
    %   - time-dependent effective rotational inertia,
    %   - time-dependent external joint torques,
    %   - desired angular position trajectory.
    %
    % INPUTS:
    %   rotInertia - a 3x3 matrix where each row contains the diagonal values of the
    %                rotational inertia for one link, e.g.,
    %                [0.02, 0.02, 0.01;
    %                 0.02, 0.02, 0.01;
    %                 0.02, 0.02, 0.01]
    %   mass       - a 1x3 vector containing the mass of each link, e.g., [1.0, 1.0, 1.0]
    %
    % OUTPUTS:
    %   time_grid    - time discretization vector (1 x N)
    %   I_eff        - time-dependent effective rotational inertia (3 x N matrix)
    %   torque_traj  - time-dependent external joint torques (3 x N matrix)
    %   angle_traj   - desired joint angular trajectory (3 x N matrix)
    %
    % Example usage:
    %   rotInertia = [0.02, 0.02, 0.01;
    %                 0.02, 0.02, 0.01;
    %                 0.02, 0.02, 0.01];
    %   mass = [1.0, 1.0, 1.0];
    %   [t, Ieff, tau, q] = PlannerFunction(rotInertia, mass);
    
        import casadi.*
        run setup.m
        %% 1. Time Discretization
        T = 2.0;    % Total duration (2 seconds)
        N = 100;    % Number of time steps
        DT = T / N;
        time_grid = linspace(0, T, N);
        
        %% 2. Create and Solve the CasADi Optimization Problem
        opti_fun = make_arm_opti();
        
        % Solve the optimization problem
        tic;
        [x, u] = opti_fun(rotInertia, mass);
        toc;
        
        X_opt = full(x);
        U_opt = full(u);
        
        NB = 3;  % Number of joints
        angle_traj = X_opt(1:NB, :);         % Joint angles [rad]
        dangle_traj = X_opt(NB+1:end, :);      % Joint velocities [rad/s]
        torque_traj = U_opt;                 % Joint torques [Nm]
        
        %% 3. Compute Time-Dependent Effective Rotational Inertia
        % Here we simply assume the effective inertia for each joint is the value from
        % the input rotational inertia along the z-axis (third element of each row).
        I_eff = zeros(NB, N);
        for i = 1:NB
            I_eff(i, :) = repmat(rotInertia(i, 3), 1, N);
        end
        
        %% 4. Plot Joint States and Inputs
        figure('Name', '3DOF Robot - States & Inputs', 'NumberTitle', 'off');
        
        subplot(1,3,1); % Joint Angles
        p1 = plot(time_grid, angle_traj(1,:), 'LineWidth', 1.5); hold on;
        p2 = plot(time_grid, angle_traj(2,:), 'LineWidth', 1.5);
        p3 = plot(time_grid, angle_traj(3,:), 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel('Angle (rad)');
        title('Joint Angles'); grid on;
        
        subplot(1,3,2); % Joint Velocities
        plot(time_grid, dangle_traj(1,:), 'LineWidth', 1.5); hold on;
        plot(time_grid, dangle_traj(2,:), 'LineWidth', 1.5);
        plot(time_grid, dangle_traj(3,:), 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel('Velocity (rad/s)');
        title('Joint Velocities'); grid on;
        
        subplot(1,3,3); % Joint Torques
        plot(time_grid, torque_traj(1,:), 'LineWidth', 1.5); hold on;
        plot(time_grid, torque_traj(2,:), 'LineWidth', 1.5);
        plot(time_grid, torque_traj(3,:), 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel('Torque (Nm)');
        title('Joint Torques'); grid on;
        
        % Add a legend outside the subplots
        hL = axes('Position',[0 0 1 1],'Visible','off');
        legend(hL, [p1, p2, p3], {'Joint1','Joint2','Joint3'}, ...
               'Orientation','horizontal','Location','northoutside');
        
        %% 5. Animate the Robot Motion
        model = create3DoFRobotModel(rotInertia, mass);
        figure('Name', '3DOF Robot Animation', 'NumberTitle','off');
        showmotion(model, time_grid, angle_traj);
        drawnow;
        
        %% 6. Compute Performance Metrics and Display Results
        avg_torque   = mean(abs(torque_traj), 2);         % Average torque [Nm]
        peak_velocity = max(abs(dangle_traj), [], 2);       % Peak velocity [rad/s]
        final_angle  = angle_traj(:, end);                  % Final joint angle [rad]
        
        results = table(avg_torque, peak_velocity, final_angle, ...
                        'VariableNames', {'Average_Torque_Nm', 'Peak_Velocity_rad_s', 'Final_Angle_rad'}, ...
                        'RowNames', {'Joint1','Joint2','Joint3'});
        disp(results);
    end
    
    %% ======== Subfunctions ========
    
    function opti_fun = make_arm_opti()
    % MAKE_ARM_OPTI Creates a CasADi function for optimizing the 3-DOF arm motion.
        import casadi.*
        
        % Time discretization
        T = 2.0;
        N = 100;
        DT = T / N;
        
        % Create the optimization problem
        opti = casadi.Opti();
        
        % Parameters: rotational inertias and masses
        inertias = opti.parameter(3,3);
        masses = opti.parameter(3);
        
        % Create the robot model
        model = create3DoFRobotModel(inertias, masses);
        NB = model.NB;
        disp(['Robot has NB = ', num2str(NB), ' joints.']);
        
        % Decision variable: x = [q; dq] where q are the joint angles and dq the velocities (size 6 x N)
        x = opti.variable(NB*2, N);
        % Control input: torque u (size 3 x N)
        u = opti.variable(NB, N);
        
        % Dynamics constraints using explicit Euler integration
        for k = 1:(N-1)
            xk = x(:, k);
            uk = u(:, k);
            qk = xk(1:NB);
            dqk = xk(NB+1:end);
            
            % Simple forward dynamics (decoupled double integrator model)
            ddqk = FDab(model, qk, dqk, uk);
            
            % Explicit Euler step
            x_next = xk + [dqk; ddqk] * DT;
            opti.subject_to( x(:, k+1) == x_next );
        end
        
        % Boundary conditions
        q_init = [0; 0; 0];
        q_goal = [pi/2; pi/4; 0];
        opti.subject_to( x(1:NB, 1) == q_init );
        opti.subject_to( x(1:NB, end) == q_goal );
        % Zero initial and final velocities
        opti.subject_to( x(NB+1:end, 1) == 0 );
        opti.subject_to( x(NB+1:end, end) == 0 );
        
        % Physical constraints: torque limits and joint angle bounds
        max_torque = 15;   % Nm
        opti.subject_to( -max_torque <= u <= max_torque );
        opti.subject_to( -pi <= x(1:NB, :) <= pi );
        
        % Objective function components:
        % (a) Control effort
        control_effort = sum(sum(u.^2));
        % (b) Smoothness (penalty on torque changes)
        smoothness_penalty = sum(sum(diff(u, 1, 2).^2));
        % (c) Final velocity penalty (optional)
        end_vel_penalty = sum(x(NB+1:end, end).^2);
        % (d) End-effector path length penalty (discrete approximation)
        end_effector_path_length = 0;
        for k = 1:(N-1)
            p_k = end_effector_position(model, x(1:NB, k));
            p_kp1 = end_effector_position(model, x(1:NB, k+1));
            end_effector_path_length = end_effector_path_length + sum1((p_kp1 - p_k).^2);
        end
        alpha_path = 5;  % weight on path length penalty
        
        % Combined objective
        J = control_effort + smoothness_penalty + end_vel_penalty + alpha_path * end_effector_path_length;
        opti.minimize(J);
        
        % Set up the solver
        opti.solver('ipopt', struct('ipopt', struct('print_level', 3)));
        
        % Convert the problem to a CasADi function
        opti_fun = opti.to_function('opti_fun', {inertias, masses}, {x, u});
    end
    
    function ddq = FDab(model, q, dq, tau)
    % FDAB Simplified forward dynamics for each joint using a decoupled double integrator model.
    % For each joint, we compute: ddq = tau ./ J_eff
    % Here, J_eff is taken as the effective rotational inertia about the z-axis.
        NB = model.NB;
        % Effective rotational inertia (extracted from the input values stored in the model)
        J_eff = [model.rotInertia_vals(1); model.rotInertia_vals(2); model.rotInertia_vals(3)];
        ddq = tau ./ J_eff;
    end
    
    function p_ee = end_effector_position(model, q_sym)
    % ENDEFFECTOR_POSITION Computes the end-effector position in the base frame.
        NB = model.NB;
        X0 = cell(NB, 1);
        
        for i = 1:NB
            [XJ, ~] = jcalc(model.jtype{i}, q_sym(i));
            if i == 1
                X0{i} = XJ * model.Xtree{i};
            else
                X0{i} = X0{i-1} * XJ * model.Xtree{i};
            end
        end
        
        X_ee = X0{NB} * model.Xtree{NB+1};
        [~, p_vec] = plux_inv(X_ee);
        p_ee = p_vec;
    end
    
    function model = create3DoFRobotModel(rotInertia_arr, mass)
    % CREATE3DOFROBOTMODEL Constructs a 3-DOF robot arm model using spatial vector algebra.
    %
    % INPUTS:
    %   rotInertia_arr - a 3x3 matrix where each row contains the diagonal values of the
    %                    rotational inertia of a link.
    %   mass           - a 1x3 vector representing the mass of each link.
    %
    % OUTPUT:
    %   model - a structure with fields:
    %           NB     : Number of joints (3)
    %           parent : Parent link indices ([0, 1, 2])
    %           jtype  : Cell array of joint types (all 'Rz')
    %           Xtree  : Cell array of fixed transforms (6x6 matrices)
    %           I      : Cell array of spatial inertia matrices (6x6 matrices)
    %
    % Additionally, the field "rotInertia_vals" stores the effective rotational inertias (about z).
        
        % Convert each row into a diagonal matrix and store in a cell array
        rotInertia = {diag(rotInertia_arr(1,:)), diag(rotInertia_arr(2,:)), diag(rotInertia_arr(3,:))};
        NB = 3;
        if length(mass) ~= NB || length(rotInertia) ~= NB
            error('mass and rotInertia must each have 3 elements!');
        end
        
        % Link lengths and twist angles (adjustable)
        L = [-0.1, 0, 0.5, 0.5];   % first three are the links; the fourth is the end-effector offset
        twist_angles = [0, pi/2, pi/2];
        
        model.NB = NB;
        model.parent = [0, 1, 2];  % Joint 1 is connected to the base; joints 2 and 3 are sequentially connected
        model.jtype = cell(NB, 1);
        model.Xtree = cell(NB+1, 1);  % NB fixed transforms plus one end-effector offset
        model.I = cell(NB, 1);
        
        model.appearance.base = {'tiles', [-1 1; -1 1; 0 0], 0.5};
        
        % Store effective rotational inertias (taken from the input, z-axis values)
        model.rotInertia_vals = zeros(NB, 1);
        
        for i = 1:NB
            model.jtype{i} = 'Rz';
            % Compute fixed transform: first rotate about y-axis by the twist angle, then translate
            R = roty(twist_angles(i));
            p = [0; 0; L(i)];
            model.Xtree{i} = plux(R, p);
            % Assume center of mass is at the midpoint of the link (along x)
            com = [L(i)/2; 0; 0];
            model.I{i} = spatialInertia(mass(i), rotInertia{i}, com);
            % Record effective rotational inertia (about z)
            model.rotInertia_vals = rotInertia_arr(:, 3);
        end
        
        % End-effector offset: translation along x by L(4)
        model.Xtree{NB+1} = plux(eye(3), [0; L(4); 0]);
        
        % (Optional) Set up appearance information for links
        model.appearance.body{1} = {{'cyl', [0, 0, 0.05; 0, 0, -0.05], 0.05}};
        model.appearance.body{2} = {{'cyl', [0, 0, 0.05; 0, 0, -0.05], 0.05}, {'cyl', [0, 0, 0; L(3), 0, 0], 0.02}};
        model.appearance.body{3} = {{'cyl', [0, 0, 0.05; 0, 0, -0.05], 0.05}, {'cyl', [0, 0, 0; 0, L(4), 0], 0.02}};
    end
    
    function I_sp = spatialInertia(m, J, p)
    % SPATIALINERTIA Computes the 6x6 spatial inertia matrix.
    %
    %   I_sp = [ J + m*(S'*S),   m*S';
    %            m*S,            m*eye(3) ]
    %
    % where S is the skew-symmetric matrix of the center of mass location p.
        S = skew(p);
        I_sp = [J + m*(S'*S), m*S'; 
                m*S,        m*eye(3)];
    end
    
    function S = skew(v)
    % SKEW Returns the 3x3 skew-symmetric matrix of vector v.
    % For v = [v1; v2; v3], the output is:
    %    [   0   -v3   v2;
    %      v3     0  -v1;
    %     -v2   v1    0 ]
        S = [    0,  -v(3), v(2);
              v(3),     0, -v(1);
             -v(2),  v(1),    0];
    end
    
    function R = roty(theta)
    % ROTY Returns the 3x3 rotation matrix for a rotation about the y-axis by angle theta.
        R = [cos(theta),  0, sin(theta);
                   0,     1,    0;
             -sin(theta), 0, cos(theta)];
    end
    
    function X = plux(R, p)
    % PLUX Constructs a 6x6 spatial transform matrix from a rotation matrix R and translation vector p.
        X = [R, zeros(3);
             skew(p)*R, R];
    end
    
    function [R, p] = plux_inv(X)
    % PLUX_INV Extracts the rotation matrix R and translation vector p from a 6x6 spatial transform matrix X.
        R = X(1:3, 1:3);
        S = X(4:6, 1:3) * R';
        p = [ S(3,2); S(1,3); S(2,1) ]; % Unskew operation
    end
    
    function [XJ, S] = jcalc(type, q)
    % JCALC Computes the joint transform XJ and motion subspace S given a joint type and joint angle q.
    % Currently implemented for the 'Rz' (revolute about z-axis) joint.
        if strcmp(type, 'Rz')
            Rz = [cos(q), -sin(q), 0;
                  sin(q),  cos(q), 0;
                  0,       0,      1];
            XJ = [Rz, zeros(3);
                  zeros(3), Rz];
            S = [0; 0; 1; 0; 0; 0];
        else
            error('jcalc for joint type %s not implemented', type);
        end
    end
    
    function showmotion(model, time_grid, Q)
    % SHOWMOTION Animates the 3-DOF robot arm given the joint trajectory Q.
        NB = model.NB;
        numSteps = length(time_grid);
        % Preallocate joint positions (3 x (NB+1) x numSteps)
        joint_positions = zeros(3, NB+1, numSteps);
        
        for k = 1:numSteps
            q = Q(:, k);
            T_total = eye(4);  % Homogeneous transformation initial value
            joint_positions(:, 1, k) = T_total(1:3, 4);  % Base origin
            for i = 1:NB
                % Compute the transformation for each joint and fixed offset
                T_joint = joint_transform(q(i), model.jtype{i});
                T_tree = spatial_to_homogeneous(model.Xtree{i});
                T_total = T_total * T_joint * T_tree;
                joint_positions(:, i+1, k) = T_total(1:3, 4);
            end
        end
        
        % Animate the robot by drawing the links frame-by-frame
        figure('Name', 'Robot Animation');
        hold on; grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
        for k = 1:numSteps
            pts = squeeze(joint_positions(:,:,k));
            plot3(pts(1,:), pts(2,:), pts(3,:), 'ko-', 'LineWidth', 1.5);
            drawnow;
            pause(0.05);
        end
    end
    
    function T = spatial_to_homogeneous(X)
    % SPATIAL_TO_HOMOGENEOUS Converts a 6x6 spatial transform into a 4x4 homogeneous transform.
        [R, p] = plux_inv(X);
        T = [R, p; 0 0 0 1];
    end
    
    function T = joint_transform(q, type)
    % JOINT_TRANSFORM Returns the 4x4 homogeneous transform corresponding to a joint rotation.
        if strcmp(type, 'Rz')
            T = [cos(q), -sin(q), 0, 0;
                 sin(q),  cos(q), 0, 0;
                 0,       0,      1, 0;
                 0,       0,      0, 1];
        else
            error('joint_transform for joint type %s not implemented', type);
        end
    end
    
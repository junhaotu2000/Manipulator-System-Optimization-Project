function opti_fun = make_arm_opti()
%MAKE_ARM_OPTI makes casadi function to 
%   Detailed explanation goes here
import casadi.*

%% 2) Time Discretization
T = 2.0;   % Total duration (2 seconds)
N = 100;   % Number of time steps
DT = T / N;
time_grid = linspace(0, T, N);

%% 3) Create CasADi Optimization Problem
opti = casadi.Opti();

inertias = opti.parameter(3,3);
masses = opti.parameter(3);
lengths = opti.parameter(3);

start_p = opti.parameter(3);
end_p = opti.parameter(3);
model = create3DoFRobotModel(inertias, masses, lengths);
NB = model.NB;  % Number of joints
disp(['Robot has NB = ', num2str(NB), ' joints.']);

% Decision variable: x = [q; dq], size (8 x N) for 3-DOF robot
x = opti.variable(NB*2, N);
% Control input: torque u, size (3 x N)
u = opti.variable(NB, N);

%% 4) Dynamics Constraints - Explicit Euler Integration
for k = 1 : (N-1)
    xk  = x(:, k);
    uk  = u(:, k);

    qk  = xk(1:NB);
    dqk = xk(NB+1:end);

    % Forward dynamics (Featherstone FDab)
    ddqk = FDab(model, qk, dqk, uk);

    % Explicit Euler step
    x_next = xk + [dqk; ddqk] * DT;

    % Apply dynamics constraint
    opti.subject_to( x(:, k+1) == x_next );
end

%% 5) Boundary Conditions
q_init = start_p; 
q_goal = end_p;

opti.subject_to( x(1:NB,1)   == q_init );
opti.subject_to( x(1:NB,end) == q_goal );

% Optional: initial/final velocity = 0
opti.subject_to( x(NB+1:end,1)   == 0 );
opti.subject_to( x(NB+1:end,end) == 0 );

%% 6) Physical Constraints
% max_torque = 15;  % Max joint torque ±15Nm
% opti.subject_to( -max_torque <= u <= max_torque );

% Joint angle limits (set according to real robot spec)
opti.subject_to( -pi <= x(1:NB,:) <= pi );

%% 7) Objective Function (Multi-objective + End-Effector Path)
% (a) Control effort
control_effort = sum(sum( u.^2 ));

% (b) Smoothness: penalty on torque changes
smoothness_penalty = sum(sum( diff(u, 1, 2).^2 ));

% (c) Final velocity penalty (optional)
end_vel_penalty = sum( x(NB+1:end, end).^2 );

% (d) End-effector path length penalty
%     Approximate by discrete sum of squared Euclidean distances
end_effector_path_length = 0;
for k = 1 : (N-1)
    p_k   = end_effector_position(model, x(1:NB, k));
    p_kp1 = end_effector_position(model, x(1:NB, k+1));
    end_effector_path_length = end_effector_path_length + sum1( (p_kp1 - p_k).^2 );
end
alpha_path = 5;  % Weight on path length penalty

% Combine all cost terms
J = control_effort + smoothness_penalty + end_vel_penalty ...
    + alpha_path * end_effector_path_length;

opti.minimize(J);

%% 8) Solver Setup and Execution
opti.solver('ipopt', struct('ipopt', struct('print_level', 3)));


opti_fun = opti.to_function('opti_fun', {inertias, masses, lengths, start_p, end_p}, {x, u});
% (Optional) Initial guess
% opti.set_initial(x, 0);
% opti.set_initial(u, 0);

end
%% =============== Compute End-Effector 3D Position ===============
function p_ee = end_effector_position(model, q_sym)
    % Returns end-effector (with offset) position in base frame
    import casadi.*
    NB_ = model.NB;
    Xup = cell(NB_,1);
    X0 = cell(NB_,1);  
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
    
    % model.Xtree{NB_+1} contains the end-effector offset
    X_ee = X0{NB_} * model.Xtree{NB_+1};

    % Extract translation vector from 6x6 transform matrix
    [R_ee, p_vec] = plux_inv(X_ee);
    p_ee = p_vec;  % 3D position in base frame
end

function [R, p] = plux_inv(X)
    % Extract rotation and translation from 6x6 spatial transform matrix
    R = X(1:3, 1:3);
    S = X(4:6, 1:3) * R'; 
    p = [ S(3,2); S(1,3); S(2,1) ];  % Unskew
end


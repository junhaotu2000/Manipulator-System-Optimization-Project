
%% 入口函数：运行力量学优化并绘图，返回最优解和最优值
function [x_opt, u_opt, Jopt] = run_arm_opti()
    %RUN_ARM_OPTI 载入参数，使用 fmincon 优化并绘图

    %% 1. 机械设计
    [rotInertia_arr, mass_out, L] = MechOptimization();
    fprintf('Mechanical Optimization Complete.\n');

    %% 2. 轨迹优化（fmincon）
    model   = create3DoFRobotModel(rotInertia_arr'*1e-6, mass_out, L/1000);
    NB      = model.NB;   N = 100;   T = 2.0;   DT = T/N;
    q_init  = [0;0;0];
    q_goal  = [pi/2;pi/4;0];
    tau_max = 15;
    alpha   = 5;

    % 初始猜测
    z0 = zeros((2*NB + NB)*N, 1);
    opts = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'Display','iter', ...
        'MaxFunctionEvaluations',1e6);

    % 调用 fmincon
    [z_opt, Jopt] = fmincon( ...
        @(z) cost_fun(z,NB,N,DT,model,alpha), ...  % 目标函数
        z0, ...                                  % 初始猜测
        [],[],[],[], ...                         % 线性约束无
        [],[], ...                                % 无显式上下界
        @(z) nonlcon(z,NB,N,DT,model,q_init,q_goal,tau_max), ... % 非线性约束
        opts );
    fprintf('Optimization Finished. J = %.4e\n', Jopt);

    % 拆解结果
    x_opt = reshape(z_opt(1:2*NB*N), 2*NB, N);
    u_opt = reshape(z_opt(2*NB*N+1:end), NB, N);

    % 绘图
    t_val = linspace(0, T, N);
    Q     = x_opt(1:NB, :);
    dQ    = x_opt(NB+1:end, :);

    figure('Name','States & Inputs','NumberTitle','off'); clf;
    subplot(1,3,1);
    plot(t_val, Q','LineWidth',1.5);
    xlabel('Time (s)'); ylabel('Angle (rad)'); title('Joint Angles'); grid on;

    subplot(1,3,2);
    plot(t_val, dQ','LineWidth',1.5);
    xlabel('Time (s)'); ylabel('Velocity (rad/s)'); title('Joint Velocities'); grid on;

    subplot(1,3,3);
    plot(t_val, u_opt','LineWidth',1.5);
    xlabel('Time (s)'); ylabel('Torque (Nm)'); title('Joint Torques'); grid on;
end

%% 目标函数
function J = cost_fun(z, NB, N, DT, model, alpha)
    x = reshape(z(1:2*NB*N), 2*NB, N);
    u = reshape(z(2*NB*N+1:end), NB, N);
    control_effort = sum(u(:).^2);
    smoothness     = sum(diff(u,1,2).^2,'all');
    end_vel        = sum(x(NB+1:end,end).^2);
    path_penalty   = 0;
    for k = 1:N-1
        p1 = ee_pos(model, x(1:NB,k));
        p2 = ee_pos(model, x(1:NB,k+1));
        path_penalty = path_penalty + sum((p2-p1).^2);
    end
    J = control_effort + smoothness + end_vel + alpha * path_penalty;
end

%% 非线性约束函数
function [c, ceq] = nonlcon(z, NB, N, DT, model, q0, qf, tau_max)
    x = reshape(z(1:2*NB*N), 2*NB, N);
    u = reshape(z(2*NB*N+1:end), NB, N);
    ceq_dyn = zeros(2*NB*(N-1),1);
    for k = 1:N-1
        qk  = x(1:NB,k); dqk = x(NB+1:end,k);
        ddq = FDab(model, qk, dqk, u(:,k));
        h   = x(:,k) + [dqk; ddq]*DT - x(:,k+1);
        ceq_dyn((k-1)*2*NB+1 : k*2*NB) = h;
    end
    ceq_bc = [ x(1:NB,1)-q0;
               x(1:NB,end)-qf;
               x(NB+1:end,1);
               x(NB+1:end,end) ];
    ceq = [ceq_dyn; ceq_bc];
    c1 = x(1:NB,:) - pi;
    c2 = -pi - x(1:NB,:);
    c3 = u - tau_max;
    c4 = -tau_max - u;
    c  = [c1(:); c2(:); c3(:); c4(:)];
end

%% 末端位姿计算
function p = ee_pos(model, q)
    NB = model.NB; X0 = cell(NB,1);
    for i = 1:NB
        [XJ,~] = jcalc(model.jtype{i}, q(i));
        if model.parent(i)==0, X0{i} = XJ*model.Xtree{i};
        else                X0{i} = X0{model.parent(i)}*XJ*model.Xtree{i}; end
    end
    X_ee = X0{NB} * model.Xtree{NB+1};
    [~, p] = plux_inv(X_ee);
end

%% 空间变换逆提取
function [R, p] = plux_inv(X)
    R = X(1:3,1:3); S = X(4:6,1:3)*R'; p = [S(3,2); S(1,3); S(2,1)];
end

%% 机器人模型构造
function model = create3DoFRobotModel(rotInertia_arr, mass, L)
    NB = 3; model.NB = NB; model.parent=[0,1,2]; types={'Rz','Rx','Rx'};
    for i=1:NB
        model.jtype{i}=types{i};
        model.Xtree{i}=plux(eye(3),[0;0;L(i)]);
        model.I{i}=mcI(mass(i),[0,0,L(i)],diag(rotInertia_arr(i,:)));
    end
    model.Xtree{NB+1}=plux(eye(3),[0;0;L(NB)]);
end

run setup.m

import casadi.*

model = autoTree(3, 1, pi/2, 1);
model.Xtree{end+1} = plux(eye(3), [1, 0,0]);

NB = model.NB;

T = 1;
N = 20;
DT = T / N;

opti = casadi.Opti();
opti.solver('ipopt')

x = opti.variable(NB * 2, N);
u = opti.variable(NB, N);


V = MX(0);
for ii = 1:N-1
    % dynamics constraint
    xi = x(:, ii);
    ui = u(:, ii);

    dxdt = dynamics(model, xi, ui);
    x(:, ii+1) = xi + dxdt * DT;

    P = ui' * ui;
    V = V + P;
end

P_init = [0; 0; 0;];
P_end = [pi/2; pi/2; 0;];

opti.subject_to(x(1:3, 1) == P_init)
opti.subject_to(x(1:3, end) == P_end)

sol = opti.solve();

%%

T_SPAN = 0:DT:T-DT;

X = sol.value(x);
U = sol.value(u);
Q = X(1:NB, :);
dQ = X(NB+1:end, :);

showmotion(model, T_SPAN, Q)

function dxdt = dynamics(model, x, u)
    nb = model.NB;
    q = x(1:nb);
    dq = x(nb+1: end);

    % [H, C] = HandC(model, q, dq);
    % H_inv = H^-1;
    % B = eye(nb);
    
    ddq = FDab(model, q, dq, u);
 
    dxdt = [dq; ddq];
end


function P = end_pos(model, x)
    nb = model.NB;
    q = x(1:nb);
    % dq = x(nb+1:end);
    
    dim = size(model.Xtree{1});
    
    if length(model.Xtree) == nb
        Xend = eye(dim);
    else
        Xend = model.Xtree{nb+1};
    end
    
    Xk = eye(dim);
    for k = 1:nb
        [Xj, ~] = jcalc(model.jtype{k}, q(k));
        % J_b(:,k) = S;
        Xup = Xj * model.Xtree{k};
        Xk = Xup * Xk;
    end
    % Xk
    Xend = Xend * Xk;
    [~, r] = plux(Xend);
    P = r;
end
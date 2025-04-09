function MechOptimization()
    clc; clear;
    % ========== 1. GA Parameter Setup ==========
    nLinks = 3;  
    % Design variable x: [D1-D3, t1-t3, L1-L3, mat1-mat3]
    % where mat_i ∈ {1,2,3,4} indicates material type
    nCont = nLinks * 3;       % Continuous variables: D, t, L
    nVar  = nCont + nLinks;   % Total number of variables
    popSize = 300;            % Population size
    nCrit   = 50;             % Number of individuals initialized with midpoints
    numRuns = 1;              % Number of GA runs

    % Variable bounds (mm; materials as integers 1–4)
    lb = [ones(1,nLinks)*20,  ones(1,nLinks)*5,  ones(1,nLinks)*50,  ones(1,nLinks)*1];
    ub = [ones(1,nLinks)*200, ones(1,nLinks)*10, ones(1,nLinks)*800, ones(1,nLinks)*4];

    % Integer indices (for material selection)
    IntCon = (nCont+1):nVar;

    % ========== 2. Initial Population Setup ==========
    initialPop = zeros(popSize, nVar);
    midpoint = (lb + ub)/2;
    for i = 1:nCrit
        initialPop(i,:) = midpoint;
    end
    for i = (nCrit+1):popSize
        initialPop(i,:) = lb + rand(1, nVar) .* (ub - lb);
        initialPop(i, IntCon) = floor(initialPop(i, IntCon));
    end

    % ========== 3. GA Options ==========
    optionsGA = optimoptions('ga', ...
        'Display','iter', ...
        'PopulationSize', popSize, ...
        'MaxGenerations', 600, ...
        'PlotFcn', {@gaplotbestf, @gaplotstopping}, ...
        'InitialPopulationMatrix', initialPop);

    % ========== 4. Run GA and Collect Results ==========
    fprintf('Running GA optimization...\n');
    [x_opt, fval_opt] = ga(@massObjective, nVar, [], [], [], [], lb, ub, @nonlcon, IntCon, optionsGA);

    % Check constraint violations
    [c_active, ~] = nonlcon(x_opt);
    activeConstraints = abs(c_active) < 1e-5; %#ok<NASGU>

    % ========== 5. Parse and Output Final Results ==========
    fprintf('\nGA Optimization Complete!\n');
    fprintf('Optimal design variables x_opt = \n');
    disp(x_opt);
    fprintf('Objective value (Total mass) = %.4f kg\n\n', fval_opt);

    % Compute final mass and link inertias
    [totalMass, linkInertias_6x6] = computeFinalInertiaAndMass(x_opt);
    fprintf('--- Final Total Mass = %.4f kg\n', totalMass);
    fprintf('--- Link 6x6 Spatial Inertia Matrices:\n');
    for i = 1:nLinks
        fprintf('Link %d inertia matrix = \n', i);
        disp(linkInertias_6x6{i});
    end
end

function f = massObjective(x)
    [mass, ~] = mechanicalModel(x);
    f = mass;
end

function [Totalmass, mechOut] = mechanicalModel(x)
    nLinks = 3;
    D = x(1:nLinks);                 % Outer diameter (mm)
    t = x(nLinks+1:2*nLinks);        % Wall thickness (mm)
    L = x(2*nLinks+1:3*nLinks);      % Length (mm)
    m_idx = round(x(3*nLinks+1:end)); % Material indices

    linkMasses = zeros(1,nLinks);
    for i = 1:nLinks
        D_m = D(i)*1e-3;
        t_m = t(i)*1e-3;
        L_m = L(i)*1e-3;
        switch m_idx(i)
            case 1, rho = 2700;
            case 2, rho = 1600;
            case 3, rho = 4500;
            case 4, rho = 7800;
            otherwise, rho = 2700;
        end
        A = (pi/4)*((D_m)^2 - (D_m - 2*t_m)^2);
        V = A * L_m;
        linkMasses(i) = rho * V;
    end

    totalMass = 10 * sum(linkMasses) + 20;
    Totalmass = sum(linkMasses);

    mechOut.D           = D;
    mechOut.t           = t;
    mechOut.L           = L;
    mechOut.materialIdx = m_idx;
    mechOut.linkMasses  = linkMasses;
end

function [c, ceq] = nonlcon(x)
    nLinks = 3;
    D = x(1:nLinks);
    t = x(nLinks+1:2*nLinks);
    L = x(2*nLinks+1:3*nLinks);
    m_idx = round(x(3*nLinks+1:end));

    S     = 4;        % Safety factor
    M_max = 500;      % Maximum bending moment (Nm)

    c = [];
    for i = 1:nLinks
        D_i = D(i);
        t_i = t(i);
        L_i = L(i);
        D_m = D_i*1e-3;
        t_m = t_i*1e-3;
        L_m = L_i*1e-3;

        I = (pi/64)*((D_m)^4 - (D_m - 2*t_m)^4);

        switch m_idx(i)
            case 1, sigma_yield = 250e6; E = 70e9;
            case 2, sigma_yield = 500e6; E = 60e9;
            case 3, sigma_yield = 900e6; E = 110e9;
            case 4, sigma_yield = 600e6; E = 210e9;
            otherwise, sigma_yield = 250e6; E = 70e9;
        end

        sigma = (M_max * (D_m / 2)) / I;
        c(end+1) = sigma - sigma_yield / S;

        c(end+1) = 2 * t_i - D_i;

        delta = (M_max * L_m^3) / (3 * E * I);
        c(end+1) = delta - 0.001 * L_m;
    end

    c(end+1) = 1000 - sum(L);
    ceq = [];
end

function [totalMass, linkInertias_6x6] = computeFinalInertiaAndMass(x)
    [totalMass, mechOut] = mechanicalModel(x);

    D = mechOut.D;
    t = mechOut.t;
    L = mechOut.L;
    linkMassArr = mechOut.linkMasses;
    nLinks = length(D);

    rotInertiaArr = zeros(3, nLinks);
    for i = 1:nLinks
        M_i  = linkMassArr(i);
        r_o  = (D(i)/2)*1e-3;
        r_i  = r_o - (t(i)*1e-3);
        L_i  = (L(i))*1e-3;

        I_x = 0.5 * M_i * (r_o^2 + r_i^2);
        I_trans = (1/12)*M_i*(3*(r_o^2 + r_i^2) + L_i^2);

        rotInertiaArr(:, i) = [I_x; I_trans; I_trans];
    end

    linkInertias_6x6 = cell(nLinks,1);
    for i = 1:nLinks
        J = diag(rotInertiaArr(:,i));
        m_i = linkMassArr(i);
        p = [L(i)*0.5e-3; 0; 0];
        S = [   0,   -p(3),  p(2);
              p(3),     0,  -p(1);
             -p(2),  p(1),    0];
        I_sp = [J + m_i*(S'*S), m_i*S';
                m_i*S,         m_i*eye(3)];
        linkInertias_6x6{i} = I_sp;
    end
end

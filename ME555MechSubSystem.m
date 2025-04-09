%% runGA.m
clear; clc;
% rng(1234);  % Uncomment for reproducibility

%% 1. Variable Definition
nLinks = 3;  
% Design variable vector: x = [D1-D3, t1-t3, L1-L3, m1-m3]
nCont = nLinks*3;       % Number of continuous variables (D, t, L)
nVar = nCont + nLinks;  % Total number of design variables (last 3 are material selections)
popSize = 300;          % Population size
nCrit = 50;             % Number of initialized individuals with midpoint values
numRuns = 5;            % Number of GA optimization runs

% Variable bounds (D, t, L in mm; m in {1,2,3,4})
lb = [ones(1, nLinks)*20, ones(1, nLinks)*5, ones(1, nLinks)*50, ones(1, nLinks)*1];
ub = [ones(1, nLinks)*200, ones(1, nLinks)*10, ones(1, nLinks)*800, ones(1, nLinks)*4];

% Integer variable indices (material indices)
IntCon = (nCont+1):nVar;

% Initial population
initialPop = zeros(popSize, nVar);
for i = 1:nCrit
    initialPop(i,:) = (lb + ub) / 2;
end
for i = (nCrit+1):popSize
    initialPop(i,:) = lb + rand(1, nVar) .* (ub - lb);
    % Ensure integer values for material selections
    initialPop(i, nCont+1:end) = floor(initialPop(i, nCont+1:end));
end

%% 2. GA Optimization Options
optionsGA = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', popSize, ...
    'MaxGenerations', 600, ...
    'PlotFcn', {@gaplotbestf, @gaplotstopping}, ...
    'InitialPopulationMatrix', initialPop);

%% 3. Run GA Optimization and Collect Results
results = cell(numRuns, 5);  % Result storage

for i = 1:numRuns
    fprintf('Starting GA optimization - Run %d...\n', i);
    
    % Use wrapped objective function
    [x_opt, fval_opt] = ga(@massObjectiveCoupled, nVar, [], [], [], [], lb, ub, @nonlcon, IntCon, optionsGA);
    
    % Evaluate constraints at optimal point
    [c_active, ~] = nonlcon(x_opt);
    active_constraints = abs(c_active) < 1e-5;
    
    % Store results
    results{i,1} = x_opt;               % Optimal design variables
    results{i,2} = fval_opt;            % Optimal objective value (mass)
    results{i,3} = initialPop(1,:);     % Reference initial point
    
    % Retrieve mechanical outputs for coupling
    [~, mechOut] = mechanicalModel(x_opt);
    fprintf('Run %d completed.\n', i);
    fprintf('Optimal Material Vector: %s\n', mat2str(mechOut.m));
    fprintf('Link Inertias: %s\n', mat2str(mechOut.linkInertias));
    fprintf('Inertia Matrix (3×3):\n');
    disp(mechOut.inertiaMatrix);
    fprintf('Link Masses: %s\n\n', mat2str(mechOut.linkMasses));
end

%% 4. Convert Results to Table
T = table(results(:,1), results(:,2), results(:,3), ...
    'VariableNames', {'x_opt', 'f_opt', 'InitialPoint'});
disp(T);

%% Encapsulated Mechanical Subsystem

% 4.1 Encapsulated mechanical model (3D coupled version)
function [mass, mechOut] = mechanicalModel(x)
    nLinks = 3;
    % Extract design variables
    D = x(1:nLinks);                           % Outer diameters (mm)
    t = x(nLinks+1:2*nLinks);                  % Wall thicknesses (mm)
    L = x(2*nLinks+1:3*nLinks);                % Link lengths (mm)
    m = round(x(3*nLinks+1:3*nLinks+nLinks));  % Material indices (1×3)
    
    % Initialize outputs
    totalMass = 0;
    linkInertias = zeros(1, nLinks);  % Sectional inertias (m^4)
    linkMasses  = zeros(1, nLinks);   % Link masses (kg)
    
    for i = 1:nLinks
        % Unit conversion: mm → m
        D_m = D(i) * 0.001;
        t_m = t(i) * 0.001;
        L_m = L(i) * 0.001;
        
        % Material density (kg/m³)
        switch m(i)
            case 1, rho = 2700;
            case 2, rho = 1600;
            case 3, rho = 4500;
            case 4, rho = 7800;
            otherwise, rho = 2700;
        end
        
        % Cross-sectional area and volume
        A = (pi/4)*((D_m)^2 - (D_m - 2*t_m)^2);
        V = A * L_m;
        mass_i = rho * V;
        totalMass = totalMass + mass_i;
        
        % Sectional moment of inertia I (m^4)
        I = (pi/64) * ((D_m)^4 - (D_m - 2*t_m)^4);
        linkInertias(i) = I;
        linkMasses(i) = mass_i;
    end
    
    % Add estimated joint and auxiliary mass
    totalMass = 10 * totalMass + 20;
    mass = totalMass;
    
    % Construct output struct for subsystem coupling
    mechOut.D = D;
    mechOut.t = t;
    mechOut.L = L;
    mechOut.m = m;
    mechOut.rho = [];  % Optionally include per-link densities
    mechOut.linkMasses = linkMasses;
    mechOut.linkInertias = linkInertias;
    mechOut.inertiaMatrix = diag(linkInertias);  % 3×3 diagonal inertia matrix
end

% 4.2 Wrapped objective function (calls mechanical model)
function f = massObjectiveCoupled(x)
    [mass, ~] = mechanicalModel(x);
    f = mass;
end

% 4.3 Nonlinear constraints (3D version)
function [c, ceq] = nonlcon(x)
    nLinks = 3;
    D = x(1:nLinks);                           % Outer diameters (mm)
    t = x(nLinks+1:2*nLinks);                  % Wall thicknesses (mm)
    L = x(2*nLinks+1:3*nLinks);                % Link lengths (mm)
    m = round(x(3*nLinks+1:3*nLinks+nLinks));  % Material indices (1×3)
    
    S = 4;           % Safety factor
    M_max = 500;     % Maximum bending moment (Nm)
    
    c = [];
    for i = 1:nLinks
        D_i = D(i);
        t_i = t(i);
        L_i = L(i);
        
        % Convert to meters
        D_m = D_i * 0.001;
        t_m = t_i * 0.001;
        L_m = L_i * 0.001;
        
        % Moment of inertia
        I = (pi/64) * ((D_m)^4 - (D_m - 2*t_m)^4);
        
        % Material properties
        switch m(i)
            case 1, sigma_yield = 250e6; E = 70e9;
            case 2, sigma_yield = 500e6; E = 60e9;
            case 3, sigma_yield = 900e6; E = 110e9;
            case 4, sigma_yield = 600e6; E = 210e9;
            otherwise, sigma_yield = 250e6; E = 70e9;
        end
        
        % Bending stress (Pa)
        sigma = (M_max * (D_m/2)) / I;
        % Stress constraint
        c(end+1) = sigma - (sigma_yield / S);
        
        % Manufacturability constraint: inner diameter > 0 → 2t - D <= 0
        c(end+1) = 2*t_i - D_i;
        
        % Deflection constraint: δ = (M_max * L³) / (3 * E * I)
        delta = (M_max * L_m^3) / (3 * E * I);
        c(end+1) = delta - 0.001 * L_m;
    end
    
    % Workspace constraint: total link length ≥ 1000 mm
    c(end+1) = 1000 - sum(L);
    
    ceq = [];
end

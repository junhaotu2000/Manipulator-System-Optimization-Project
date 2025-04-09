function [rotInertia_arr, mass_out] = MechOptimization()
% runOptimizedMechanicalArm integrates the GA optimization and inertia calculation for the mechanical subsystem.
%
% Outputs:
%   rotInertia_arr - 3x3 matrix, where each column corresponds to a link's [I_x; I_y; I_z]
%   mass_out       - 1x3 vector of actual link masses
%
% Overall workflow:
%   1. Use a genetic algorithm to optimize the design vector x. The design vector includes:
%          D: link outer diameter (mm)
%          t: link wall thickness (mm)
%          L: link length (mm)
%          m: material selection (integer 1 to 4, corresponding to different densities)
%   2. In the function mechanicalModel, compute each link's actual mass (linkMasses) and inertia (linkInertias) from x.
%   3. Use convertLinkParams (internally encapsulated within the core of create3DoFRobotModel) to compute
%      each link's rotational inertia parameters (rotInertia_arr) and the actual mass vector (mass_out) based on D, t, L, and linkMasses.
%
% During execution, the GA optimization process and some mechanical subsystem information after optimization are displayed in the command window,
% and the two required outputs are returned at the end.

    %% 1. GA Optimization Settings
    nLinks = 3;  
    % Design vector x = [D1, D2, D3, t1, t2, t3, L1, L2, L3, m1, m2, m3]
    nCont = nLinks * 3;       % Number of continuous variables (D, t, L)
    nVar = nCont + nLinks;    % Total number of design variables (last 3 are for material selection)
    popSize = 300;            % Population size
    nCrit = 50;               % Number of individuals in the initial population set to the midpoint value
    numRuns = 3;              % Number of GA optimization runs; adjust as needed

    % Variable bounds (units: mm for D, t, L; m is an index in the range 1 to 4)
    lb = [ones(1, nLinks)*20, ones(1, nLinks)*5, ones(1, nLinks)*50, ones(1, nLinks)*1];
    ub = [ones(1, nLinks)*200, ones(1, nLinks)*10, ones(1, nLinks)*800, ones(1, nLinks)*4];
    % Integer variable indices
    IntCon = (nCont+1):nVar;

    % Build the initial population
    initialPop = zeros(popSize, nVar);
    for i = 1:nCrit
        initialPop(i,:) = (lb + ub) / 2;
    end
    for i = (nCrit+1):popSize
        initialPop(i,:) = lb + rand(1, nVar) .* (ub - lb);
        % Material selection must be an integer
        initialPop(i, nCont+1:end) = floor(initialPop(i, nCont+1:end));
    end

    % GA optimization options
    optionsGA = optimoptions('ga', ...
        'Display', 'iter', ...
        'PopulationSize', popSize, ...
        'MaxGenerations', 600, ...
        'PlotFcn', {@gaplotbestf, @gaplotstopping}, ...
        'InitialPopulationMatrix', initialPop);

    %% 2. Running GA Optimization
    results = cell(numRuns, 2);  % Store the optimal design vector and objective values
    best_x = [];
    best_fval = inf;
    for i = 1:numRuns
        fprintf('Starting GA optimization - Run %d...\n', i);
        [x_opt, fval_opt] = ga(@massObjectiveCoupled, nVar, [], [], [], [], lb, ub, @nonlcon, IntCon, optionsGA);
        
        % Optional: Check constraint function and other details here
        
        % Update the global optimal solution
        if fval_opt < best_fval
            best_fval = fval_opt;
            best_x = x_opt;
        end
        
        fprintf('Run %d completed. Optimal objective (total mass): %.4f\n', i, fval_opt);
        
        % Call the mechanical subsystem model to obtain output parameters (including D, t, L, and each link's actual mass)
        [~, mechOut] = mechanicalModel(x_opt);
        fprintf('Optimal Material Vector (material selection indices): %s\n', mat2str(mechOut.m));
        fprintf('Individual link masses (kg): %s\n\n', mat2str(mechOut.linkMasses));
        fprintf('Individual link lengths (mm): %s\n\n', mat2str(mechOut.L));
        results{i,1} = x_opt;
        results{i,2} = fval_opt;
    end

    %% 3. Calculate Rotational Inertia Matrix and Link Mass Vector
    % Retrieve the optimal design parameters from the mechanical subsystem model
    [~, mechOut] = mechanicalModel(best_x);
    % Note: convertLinkParams expects D, t, L in mm and mass in kg.
    %       Here we use the actual link masses computed by the mechanical subsystem.
    [rotInertia_arr, mass_out] = convertLinkParams(mechOut.D, mechOut.t, mechOut.L, mechOut.linkMasses);
    
    %% Display Final Outputs
    fprintf('\nFinal Outputs:\n');
    disp('rotInertia_arr (each column is [I_x; I_y; I_z] in kg*mm^2):');
    disp(rotInertia_arr);
    disp('mass_out (link masses in kg):');
    disp(mass_out);
    fprintf('Individual link lengths (mm): %s\n\n', mat2str(mechOut.L));

end

%% ----------------- Subfunctions Section ------------------

function f = massObjectiveCoupled(x)
    % massObjectiveCoupled wraps the objective function and returns the total mass computed by the mechanical subsystem.
    [mass, ~] = mechanicalModel(x);
    f = mass;
end

function [mass, mechOut] = mechanicalModel(x)
    % mechanicalModel calculates the output parameters of the mechanical subsystem based on the design vector x.
    % Input:
    %   x - Design vector [D, t, L, m],
    %       where D: outer diameter (mm), t: wall thickness (mm), L: link length (mm), 
    %       and m: material selection (integer 1 to 4)
    %
    % Outputs:
    %   mass   - Total mass (including estimated masses for joints and auxiliary parts) used for objective feedback
    %   mechOut- Structure containing optimized design parameters and computed link parameters (D, t, L, m, linkMasses, linkInertias, inertiaMatrix)
    
    nLinks = 3;
    % Extract design variables
    D = x(1:nLinks);                           % Outer diameter (mm)
    t = x(nLinks+1:2*nLinks);                    % Wall thickness (mm)
    L = x(2*nLinks+1:3*nLinks);                  % Link length (mm)
    m_index = round(x(3*nLinks+1:3*nLinks+nLinks));  % Material selection (integer, 1x3)
    
    totalMass = 0;
    linkInertias = zeros(1, nLinks);   % Cross-sectional inertias for further estimation (in m^4), note that these differ from later inertia formulas
    linkMasses  = zeros(1, nLinks);    % Actual mass of each link (kg)
    
    for i = 1:nLinks
        % Unit conversion: mm → m
        D_m = D(i) * 0.001;
        t_m = t(i) * 0.001;
        L_m = L(i) * 0.001;
        
        % Determine density (kg/m³) based on material selection
        switch m_index(i)
            case 1, rho = 2700;
            case 2, rho = 1600;
            case 3, rho = 4500;
            case 4, rho = 7800;
            otherwise, rho = 2700;
        end
        
        % Calculate cross-sectional area and volume (hollow circular tube)
        A = (pi/4)*((D_m)^2 - (D_m - 2*t_m)^2);
        V = A * L_m;
        mass_i = rho * V;
        totalMass = totalMass + mass_i;
        
        % Calculate cross-sectional inertia (for reference in structural load estimation)
        I = (pi/64) * ((D_m)^4 - (D_m - 2*t_m)^4);
        linkInertias(i) = I;
        linkMasses(i) = mass_i;
    end
    
    % Include estimated mass for joints and auxiliary components (scaling factor and offset)
    totalMass = 10 * totalMass + 20;
    mass = totalMass;
    
    % Construct the mechanical subsystem output structure
    mechOut.D = D;
    mechOut.t = t;
    mechOut.L = L;
    mechOut.m = m_index;          % Material selection (original design variable)
    mechOut.linkMasses = linkMasses;
    mechOut.linkInertias = linkInertias;
    mechOut.inertiaMatrix = diag(linkInertias);  % Create a 3x3 diagonal inertia matrix (for reference)
end

function [c, ceq] = nonlcon(x)
    % nonlcon defines the nonlinear constraints to ensure the design meets mechanical performance and manufacturing requirements.
    nLinks = 3;
    D = x(1:nLinks);                           % Outer diameter (mm)
    t = x(nLinks+1:2*nLinks);                    % Wall thickness (mm)
    L = x(2*nLinks+1:3*nLinks);                  % Link length (mm)
    m_index = round(x(3*nLinks+1:3*nLinks+nLinks));  % Material selection
    
    S = 4;           % Safety factor
    M_max = 500;     % Maximum bending moment (Nm)
    
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
        I = (pi/64) * ((D_m)^4 - (D_m - 2*t_m)^4);
        
        % Material properties
        switch m_index(i)
            case 1, sigma_yield = 250e6; E = 70e9;
            case 2, sigma_yield = 500e6; E = 60e9;
            case 3, sigma_yield = 900e6; E = 110e9;
            case 4, sigma_yield = 600e6; E = 210e9;
            otherwise, sigma_yield = 250e6; E = 70e9;
        end
        
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

function [rotInertia_arr, mass_out] = convertLinkParams(D, t, L, m)
    % convertLinkParams computes the rotational inertia and mass for each link based on cylindrical parameters.
    %
    % Inputs:
    %   D - 1x3 vector of outer diameters (mm)
    %   t - 1x3 vector of wall thicknesses (mm)
    %   L - 1x3 vector of link lengths (mm)
    %   m - 1x3 vector of link masses (kg)
    %
    % Outputs:
    %   rotInertia_arr - 3x3 matrix, each column is [I_x; I_y; I_z], where:
    %         I_x = 0.5*m*(r_outer^2 + r_inner^2)    (about the link axis)
    %         I_y = I_z = (1/12)*m*(3*(r_outer^2 + r_inner^2) + L^2)
    %   mass_out       - 1x3 mass vector (same as input m)

    if numel(D) ~= 3 || numel(t) ~= 3 || numel(L) ~= 3 || numel(m) ~= 3
        error('All input vectors D, t, L, and m must have exactly 3 elements!');
    end

    rotInertia_arr = zeros(3, 3);
    mass_out = m;

    for i = 1:3
        r_outer = D(i) / 2;
        r_inner = r_outer - t(i);
        if r_inner < 0
            error('The wall thickness for link %d is too large, resulting in a negative inner radius!', i);
        end

        I_x = 0.5 * m(i) * (r_outer^2 + r_inner^2);
        I_trans = (1/12) * m(i) * (3*(r_outer^2 + r_inner^2) + L(i)^2);
        rotInertia_arr(:, i) = [I_x; I_trans; I_trans];
    end
end

% Auxiliary functions (for spatial inertia, transformation matrices, etc.) for future extensions

function I_sp = spatialInertia(m, J, p)
    % spatialInertia computes the 6x6 spatial inertia matrix.
    % I_sp = [ J + m*(S'*S),   m*S';
    %          m*S,            m*I3 ]
    S = skew(p);
    I_sp = [J + m*(S'*S), m*S';
            m*S,         m*eye(3)];
end

function S = skew(v)
    % skew returns the 3x3 skew-symmetric matrix of vector v.
    S = [    0,  -v(3),  v(2);
          v(3),     0, -v(1);
         -v(2),  v(1),    0];
end

function R = roty(theta)
    % roty returns the 3x3 rotation matrix for a rotation of theta (radians) about the y-axis.
    R = [cos(theta),  0, sin(theta);
              0,      1,      0;
         -sin(theta), 0, cos(theta)];
end

function X = plux(R, p)
    % plux constructs a 6x6 spatial transformation matrix from the rotation matrix R and translation vector p.
    X = [R, zeros(3); 
         skew(p)*R, R];
end
